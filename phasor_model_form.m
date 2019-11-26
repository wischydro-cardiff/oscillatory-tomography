function [A_div,A_omegamod,b_bccoeff] = phasor_model_form(K,Ss,xb,yb,zb,bdry_types,varargin)

% phasor_model_form: Low-level function which assembles the linear algebra
% problem for the steady-periodic groundwater flow equation in 3D on a
% finite-difference-type grid with a rectangular prism as boundary.
% Boundary conditions of constant head or no flux may be applied at each of
% the 6 surrounding boundaries, but must be constant across each face of
% the rectangular prism.
%
% [A_div,A_omegamod,b_bccoeff] = phasor_model_form(K,Ss,xb,yb,zb,bdry_types)
%   OUTPUTS:
%       A_div - (num_cells x num_cells) "Stiffness" matrix for the system,
%       physically representing flows into each given cell for a particular
%       set of boundary types
%       A_omegamod - (num_cells x num_cells) Diagonal matrix used to modify
%       the system for a particular oscillation frequency (discussed in
%       further detail below)
%       b_bccoeff - (num_cells x 6) Matrix used to calculate boundary
%       contribution to b matrix, when given head values at each of the 6
%       boundaries (low x, high x, low y, high y, low z, high z)
%
% To use these variables to solve a problem with a given source term
% operating at a given angular frequency (omega), one solves the equation:
%
%              (A_steady + A_omegamod*omega)*Phi = (inflow_vec +
%              b_bccoeff*bdry_vals)
%
% where inflow_vec is the external phasor flux into each cell (positive
% inflow_vec = source) as reshaped into a flow vector with standard
% meshgrid mapping, bdry_vals is a 6 x 1 vector of the boundary head values
% (note, if any boundary is no-flux, the corresponding values will not be
% used), and omega is the angular frequency of the test
%
% ***Notes:***
%   1) The unknowns in the 3D grid are mapped to a vector using the
% standard method returned by meshgrid, i.e. phi =
% reshape(phi_grid,num_cells,1), such that phi starts in the lower
% South-west corner of the grid and increments through y first, then
% through x, then through z.
%   2) This code can also be used to solve 2D problems by setting the 3rd
% dimension to have one layer with no flux boundaries on either side.
%   3) This code uses a block-centered setup where inter-block conductances
% are calculated assuming K_1 up to a cell boundary and K_2 after that cell
% boundary, similar to MODFLOW. This results in inter-block conductances
% that are a function of geometric averages of neighboring blocks' K values.
%
% Syntax: [A_steady,A_omegamod,b_bccoeff] =
% phasor_model_form(K,Ss,xb,yb,zb,bdry_types)
% where:
%   -[A_steady,A_omegamod,b_bccoeff] are the matrices used in the problem
%   solution, as discussed above
%   -K is the matrix of hydraulic conductivity values in the cells of the
%   model, a (num_y x num_x x num_z) matrix, as would be produced from
%   meshgridding the cell centers
%   -Ss is the matrix of Specific storage values in the cells of the model,
%   a (num_y x num_x x num_z) matrix, as would be produced from
%   meshgridding the cell centers.
%   -xb is a ((num_x + 1) x 1) monotonically increasing
%   vector containing the positions of cell boundaries along the x
%   dimension.
%   -yb and zb are similar vectors to xb, but for the y and z dimensions,
%   respectively
%   -bdry_types is a 6 x 1 vector that defines the boundary condition
%   types. If an element of bdry_types is 1, it means that that boundary is
%   "constant head" type, if it is 0, it is treated as a "no flux" boundary.
%   bdry_types = [ bdry at smallest x;
%                  bdry at largest x;
%                  bdry at smallest y;
%                  bdry at largest y;
%                  bdry at smallest z; 
%                  bdry at largest z; ]
%
% The problem dimension to be solved is num_cells, where num_cells = num_x
% * num_y * num_z. Note num_x, num_y, and num_z are the number of CELLS in
% each dimension, which is one less than the number of boundaries in each
% dimension.
%
% Code based off of a similar setup for steady-state pumping tests (Michael
% Cardiff 7/2010) and modified to extend and make more efficient for
% steady-periodic pumping tests, along with additional cleanup (Michael
% Cardiff, 6/2014, 1/2015, 1/2016)
% 

%IN PROGRESS / BETA: Option for Sy linearized water table boundary
%condition
num_reqin = 6;
Sy = [];
if nargin > num_reqin
    Sy = varargin{1};
end

num_x = numel(xb) - 1;
num_y = numel(yb) - 1;
num_z = numel(zb) - 1;
num_cells = num_x*num_y*num_z;

%Check sizes of all parameters
if (size(K,1) ~= num_y) || (size(K,2) ~= num_x) || (size(K,3) ~= num_z)
    error(['Ks must be num_y by num_x by num_z, i.e. matching the ', ...
        'number of cells in each dimension']);
end
if (size(Ss,1) ~= num_y) || (size(Ss,2) ~= num_x) || (size(Ss,3) ~= num_z)
    error(['Ss must be num_y by num_x by num_z, i.e. matching the ', ...
        'number of cells in each dimension']);
end
if ~isempty(Sy)
    if (size(Sy,1) ~= num_y) || (size(Sy,2) ~= num_x)
        error(['Sy must be a num_y by num_x matrix, i.e. matching the number of ', ...
            'cells at the top of the model'])
    end
end
if (size(bdry_types) ~= [6 1])
    error('bdry_types must be a 6 x 1 matrix');
end

xtypes = bdry_types([1 2]);
ytypes = bdry_types([3 4]);
ztypes = bdry_types([5 6]);

if (num_x == 1) && (sum(xtypes == 0) ~= 2)
    error(['If there is only one cell along a dimension, the corresponding', ...
        'boundary conditions types must be no-flow (0)']);
end
if (num_y == 1) && (sum(ytypes == 0) ~= 2)
    error(['If there is only one cell along a dimension, the corresponding', ...
        'boundary conditions types must be no-flow (0)']);
end
if (num_z == 1) && (sum(ztypes == 0) ~= 2)
    error(['If there is only one cell along a dimension, the corresponding', ...
        'boundary conditions types must be no-flow (0)']);
end

if (sum(bdry_types(1:5) == 2) > 0)
    error(['Only boundary 6 (the top in the z direction) can be a', ...
        'water table boundary']);
end

if (ztypes(2) == 2) && isempty(Sy)
    error(['A water table boundary condition has been specified, but ', ...
        'no Sy input was given'])
end

dx = xb(2:end) - xb(1:(end-1));
dy = yb(2:end) - yb(1:(end-1));
dz = zb(2:end) - zb(1:(end-1));

%Setup so that all cell size vectors are column vectors. Needs to be
%standardized to match dimensions of vectorization later.
if size(dx,1) == 1
    dx = dx';
end
if size(dy,1) == 1
    dy = dy';
end
if size(dz,1) == 1
    dz = dz';
end

%Note for all reshapes below: Due to the way reshape operates, vector of
%unknowns will start from bottom, lower-SW corner of the domain and increment
%through Y first, then X, then Z. In other words, starting from lower-SW
%corner, we travel north, then east, then up during incrementing.

%SECTION: Calculate all coefficients associated with inter-cell mass
%transfer (cell-to-cell flows)
A_coldiags = zeros(num_cells,7);

%Using cell i,j,k as the center of the stencil, derive coefficients for
%flux in from lower side of stencil (i.e., -k direction).
%Formula: flux into (i,j,k) from bottom
% = Conductance(i,j,k-.5)*(h(i,j,k-1) - h(i,j,k))
A_bottom = zeros(num_y,num_x,num_z);
for k = 2:1:(num_z)
    for j = 1:1:(num_x)
        A_bottom(:,j,k) = 2*(K(:,j,k).*K(:,j,k-1).*dy*dx(j)) ...
            ./(K(:,j,k)*dz(k-1) + K(:,j,k-1)*dz(k));
        %A_bottom(:,j,k) = 5; %For testing ordering
    end
end
A_coldiags(:,1) = reshape(A_bottom,num_cells,1);
% clear A_bottom

%Using cell i,j,k as the center of the stencil, derive coefficients for
%flux in from western side of stencil (i.e., -j direction).
%Formula: flux into (i,j,k) from west 
% = Conductance(i,j-.5,k)*(h(i,j-1,k) - h(i,j,k))
A_west = zeros(num_y,num_x,num_z);
for j = 2:1:(num_x)
    for k = 1:1:(num_z)
        A_west(:,j,k) = 2*(K(:,j,k).*K(:,j-1,k).*dy*dz(k)) ...
            ./(K(:,j,k)*dx(j-1) + K(:,j-1,k)*dx(j));
        %A_west(:,j,k) = 1; %For testing ordering
    end
end
A_coldiags(:,2) = reshape(A_west,num_cells,1);
% clear A_west

%Using cell i,j,k as the center of the stencil, derive coefficients for
%flux in from southern side of stencil (i.e., -i direction).
%Formula: flux into (i,j,k) from south
% = Conductance(i-.5,j,k)*(h(i-1,j,k) - h(i,j,k))
A_south = zeros(num_y,num_x,num_z);
for i = 2:1:(num_y)
    for k = 1:1:(num_z)
        A_south(i,:,k) = 2*(K(i,:,k).*K(i-1,:,k).*dx'*dz(k)) ...
            ./(K(i,:,k)*dy(i-1) + K(i-1,:,k)*dy(i));
        %A_south(i,:,k) = 3; %For testing ordering
    end
end
A_coldiags(:,3) = reshape(A_south,num_cells,1);
% clear A_south

%Using cell i,j,k as the center of the stencil, derive coefficients for
%flux in from northern side of stencil (i.e., +i direction).
%Formula: flux into (i,j,k) from north
% = Conductance(i+.5,j,k)*(h(i+1,j,k) - h(i,j,k))
A_north = zeros(num_y,num_x,num_z);
for i = 1:1:(num_y-1)
    for k = 1:1:(num_z)
        A_north(i,:,k) = 2*(K(i,:,k).*K(i+1,:,k).*dx'*dz(k)) ...
            ./(K(i,:,k)*dy(i+1) + K(i+1,:,k)*dy(i));
        %A_north(i,:,k) = 4; %For testing ordering
    end
end
A_coldiags(:,5) = reshape(A_north,num_cells,1);
% clear A_north

%Using cell i,j,k as the center of the stencil, derive coefficients for
%flux in from eastern side of stencil (i.e., +j direction).
%Formula: flux into (i,j,k) from east 
% = Conductance(i,j+.5,k)*(h(i,j+1,k) - h(i,j,k))
A_east = zeros(num_y,num_x,num_z);
for j = 1:1:(num_x-1)
    for k = 1:1:(num_z)
        A_east(:,j,k) = 2*(K(:,j,k).*K(:,j+1,k).*dy*dz(k)) ...
            ./(K(:,j,k)*dx(j+1) + K(:,j+1,k)*dx(j));
        %A_east(:,j,k) = 2; %For testing ordering
    end
end
A_coldiags(:,6) = reshape(A_east,num_cells,1);
% clear A_east

%Using cell i,j,k as the center of the stencil, derive coefficients for
%flux in from top side of stencil (i.e., +k direction).
%Formula: flux into (i,j,k) from bottom
% = Conductance(i,j,k+.5)*(h(i,j,k+1) - h(i,j,k))
A_top = zeros(num_y,num_x,num_z);
for k = 1:1:(num_z-1)
    for j = 1:1:(num_x)
        A_top(:,j,k) = 2*(K(:,j,k).*K(:,j,k+1).*dy*dx(j)) ...
            ./(K(:,j,k)*dz(k+1) + K(:,j,k+1)*dz(k));
        %A_top(:,j,k) = 6; %For testing ordering
    end
end
A_coldiags(:,7) = reshape(A_top,num_cells,1);
% clear A_top

%SECTION: Calculate all coefficients associated with constant head boundary
%conditions. If only one cell is included in any dimension, then the
%corresponding boundary conditions must be no-flux.
b_bccoeff = sparse(num_cells,6);
if num_x > 1
    bc_coeff_west = zeros(num_y,num_x,num_z);
    if xtypes(1) == 1
        for k = 1:1:num_z
            bc_coeff_west(:,1,k) = 2.*K(:,1,k).*dy*dz(k)/dx(1);
        end
    end
    b_bccoeff(:,1) = reshape(bc_coeff_west,num_cells,1);
    % clear bc_coeff_west
    bc_coeff_east = zeros(num_y,num_x,num_z);
    if xtypes(2) == 1
        for k = 1:1:num_z
            bc_coeff_east(:,end,k) = 2.*K(:,end,k).*dy*dz(k)/dx(end);
        end
    end
    b_bccoeff(:,2) = reshape(bc_coeff_east,num_cells,1);
    % clear bc_coeff_east
end

if num_y > 1
    bc_coeff_south = zeros(num_y,num_x,num_z);
    if ytypes(1) == 1
        for k = 1:1:num_z
            bc_coeff_south(1,:,k) = 2.*K(1,:,k).*dx'*dz(k)/dy(1);
        end
    end
    b_bccoeff(:,3) = reshape(bc_coeff_south,num_cells,1);
    % clear bc_coeff_south
    bc_coeff_north = zeros(num_y,num_x,num_z);
    if ytypes(2) == 1
        for k = 1:1:num_z
            bc_coeff_north(end,:,k) = 2.*K(end,:,k).*dx'*dz(k)/dy(end);
        end
    end
    b_bccoeff(:,4) = reshape(bc_coeff_north,num_cells,1);
    % clear bc_coeff_north
end

if num_z > 1
    bc_coeff_bottom = zeros(num_y,num_x,num_z);
    if ztypes(1) == 1
        for j = 1:1:num_x
            bc_coeff_bottom(:,j,1) = 2.*K(:,j,1).*dy*dx(j)/dz(1);
        end
    end
    b_bccoeff(:,5) = reshape(bc_coeff_bottom,num_cells,1);
    % clear bc_coeff_bottom
    bc_coeff_top = zeros(num_y,num_x,num_z);
    if ztypes(2) == 1
        for j = 1:1:num_x
            bc_coeff_top(:,j,end) = 2.*K(:,j,end).*dy*dx(j)/dz(end);
        end
    end
    b_bccoeff(:,6) = reshape(bc_coeff_top,num_cells,1);
    % clear bc_coeff_top
end

%Coefficient terms for the A main diagonal comes from negative sum of all
%other diagonals (fluxes from surrounding active cells) along with fluxes
%from constant head boundary conditions. After this step, A is properly
%defined in terms of its diagonals.
%Main diagonal is negative sum of all coefficients, before considering
%boundary conditions

%Based on coefficients calculated above, create a matrix where each column
%represents a given diagonal of the A matrix.
A_coldiags(:,4) = sum(A_coldiags,2) + b_bccoeff*ones(6,1); 
A_coldiags(:,[1:3,5:7]) = -A_coldiags(:,[1:3,5:7]);
%Ordering note for A_coldiags:
%[bottom west south center north east top]

%Create A_steady based on the diagonals computed above.
A_div = spalloc(num_cells,num_cells,7*num_cells);
%spdiags does diagonals in a somewhat strange way. This results in
%diagonals where the first elements of the input vector are cut off if
%placed above the main diagonal, and the last elements of the input vector
%are cut off if placed below the main daigonal. In order to get the
%behavior we want (which is exactly opposite), diagonals are placed
%opposite to where they should be, and then transposed.
A_div = spdiags(A_coldiags,...
    [num_y*num_x num_y 1  0 -1 -num_y -num_y*num_x],A_div);
A_div = A_div.';
%Immediately southern cells actually end up on first sub-diagonal
%Immediately northern cells actually end up on first super-diagonal, etc.

%Section: Modifications to A dependent on frequency
A_omegamod_coldiags = zeros(num_cells,2);
%Ordering note for A_omegamod_coldiags:
%[bottom center]

%Contribution from specific storage to main daigonal
[dxg, dyg, dzg] = meshgrid(dx,dy,dz);
dxdydz = dxg.*dyg.*dzg;
clear dxg dyg dzg
A_omegamod_center_Ss = 1i.*Ss.*dxdydz;
clear dxdydz

%Terms in A_omegamod matrix in the case of a linearized water table boundary.
%NOTE***: Not supported. Use of specific yield (linearized water table
%boundary condition) should be considered as "in beta"
A_omegamod_center_lw = zeros(num_y,num_x,num_z);
A_omegamod_bottom_lw = zeros(num_y,num_x,num_z);
if ztypes(2) == 2
    for i = 1:1:num_y
        for j = 1:1:num_x
            A_omegamod_center_lw(i,j,end) = -1i*Sy(i,j)*dx(j)*dy(i)*(1+(dz(end)/(dz(end) + dz(end-1))));
            A_omegamod_bottom_lw(i,j,end-1) = 1i*Sy(i,j)*dx(j)*dy(i)*(dz(end)/(dz(end) + dz(end-1)));
        end
    end
end
A_omegamod_coldiags(:,1) = reshape(A_omegamod_bottom_lw,num_cells,1);
A_omegamod_coldiags(:,2) = reshape(A_omegamod_center_Ss + A_omegamod_center_lw,num_cells,1);

%Create A_steady based on the diagonals computed above.
A_omegamod = spalloc(num_cells,num_cells,2*num_cells);
%Again, spdiags does diagonals in a somewhat strange way. This results in
%diagonals where the first elements of the input vector are cut off if
%placed above the main diagonal, and the last elements of the input vector
%are cut off if placed below the main daigonal. In order to get the
%behavior we want (which is exactly opposite), diagonals are placed
%opposite to where they should be, and then transposed.
A_omegamod = spdiags(A_omegamod_coldiags,...
    [num_y*num_x 0],A_omegamod);
%NOTE: Here .' **MUST** be used instead of ' In MATLAB the symbol ' gives
%the conjugate transpose (swapping sign of imaginary part), which is not
%what we want. This was a difficult to catch bug!
A_omegamod = A_omegamod.';
%Immediately southern cells actually end up on first sub-diagonal
%Immediately northern cells actually end up on first super-diagonal, etc.
