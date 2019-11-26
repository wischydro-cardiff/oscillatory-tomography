function [sim_obs,varargout] = phasor_model_obssens...
    (omega,obs_list,obsweight_list,inflows_list,bdry_vals,K,Ss,...
    xb,yb,zb,bdry_types,varargin)

% phasor_model_obssens: Low-level function which sets up a numerical model
% for a steady-periodic problem and produces simulated observations, along
% with the full field variable solution and adjoint sensitivities, if
% requested.
%
% Notes: 
% 1) The program is built to calculate observations, field variables, and
% sensitivities for a single testing frequency, though observations of
% multiple stimulation scenarios (pumping tests) can be evaluated at the
% same time.
% 2) The program likewise assumes that the same boundary condition types
% and values apply for all stimulation scenarios.
% 3) An observation consists of a particular "obsweight" (a weighting inner
% product on the Phi field variable) applied to the results of Phi under a
% given stimulation ("fluxes") scenario. 
%
% Syntax: [sim_obs, [Phi], [H]] =
% phasor_model_obssens(omega, obs_list, obsweight_list,
% inflows_list, bdry_vals, K, Ss, xb, yb, zb, bdry_types,[mat_solver])
% where:
%   OUTPUTS:
%      -sim_obs is a vector of simulated observations (num_obs x 1). Note this
%      will generally include observations performed during tests with
%      different stimulation locations (though all at the same angular
%      frequency)
%      -Phi is the solution for the field variable across each stimulation
%      scenario, having dimensions (num_cells x num_stims)
%      -H is the sensitivity of each observation to the spatial distribution
%      of aquifer parameters. The size of H is (num_obs x num_cells*2), where
%      the first num_cells columns of H represent sensitivity to K and the
%      second num_cells columns represent sensitivity to Ss.
%   INPUTS:
%      -omega (scalar) is the angular frequency of stimulation used for all
%      steady-periodic stimulation (pumping) scenarios 
%      -obs_list (num_obs x 2) defines each observation in terms of the
%      stimulation scenario under which the observation takes place (column 1)
%      and the obsweights that are applied to obtain the observation (column
%      2)
%      -obsweight_list (num_cells x num_obstypes) is a matrix where each
%      column defines how a solution for Phi is converted into a particular
%      type of measurement. For example, for a point measurement right at a
%      grid cell center, obsweight_list would have a 1 exactly at the grid
%      cell center and a 0 everywhere else, such that the observation is given
%      by: obsweight_list(:,obstype)'*Phi(:,stim_type)
%      -inflows_list (num_cells x num_stims) is a matrix where each column
%      defines the flows into the domain for a given stimulation scenario.
%      -bdry_vals (6 x 1) contains the boundary values associated with each of
%      the boundaries (same format as bdry_types). Note that these values only
%      end up being used if the associated bdry_type = 1, representing a
%      specified value boundary.
%      -All other inputs (bdry_vals, K, Ss, xb, yb, zb, bdry_types) are as
%      defined in phasor_model_form
%      -mat_solver (optional) is an anonymous function capable of solving
%      for Phi in the numerical problem A*Phi = b, where A contains complex
%      numbers (num_cells x num_cells) and b is a matrix (num_cells *
%      num_stims)
%
% Code by Michael Cardiff, 2015-2016

num_reqin = 11;

%12th input, if supplied, is function used to invert key matrices. Matlab \
%operator is the default.
matinv_func = 'mldivide';
if nargin > (num_reqin)
    matinv_func = varargin{2};
end

%13th input, if supplied, is treated as Sy matrix. NOT YET supported
Sy = [];
if nargin > (num_reqin + 1)
    Sy = varargin{1};
end


if isempty(Sy)
    [A_steady,A_omegamod,b_bccoeff] = phasor_model_form(K,Ss,xb,yb,zb,bdry_types);
else
    [A_steady,A_omegamod,b_bccoeff] = phasor_model_form(K,Ss,xb,yb,zb,bdry_types,Sy);
end
%Form A for the particular testing frequency being carried out. Any errors
%should be caught by the steadyperiodic_poisson_form code
A = (A_steady + A_omegamod*omega);


num_obs = size(obs_list,1);
num_x = numel(xb) - 1;
num_y = numel(yb) - 1;
num_z = numel(zb) - 1;
num_cells = num_x*num_y*num_z;
num_stims = size(inflows_list,2);

xtypes = bdry_types([1 2]);
ytypes = bdry_types([3 4]);
ztypes = bdry_types([5 6]);

dx = xb(2:end) - xb(1:(end-1));
dy = yb(2:end) - yb(1:(end-1));
dz = zb(2:end) - zb(1:(end-1));

%Setup so that all cell size vectors are column vectors. Needs to be
%standardized to match dimensions of vectorization later.
if size(dx,1) == 1
    dx = dx.';
end
if size(dy,1) == 1
    dy = dy.';
end
if size(dz,1) == 1
    dz = dz.';
end

%Check to see if boundary condition values are given for each stimulation
%setup. If not, then just assume that it is the same for every stimulation
%setup.
if size(bdry_vals,2) == 1
    bdry_vals = repmat(bdry_vals,1,num_stims);
elseif size(bdry_vals,2) ~= num_stims
    error('bdry_vals must either by 6 x 1 or 6 x num_stims');
end

if size(bdry_vals,1) ~= 6
    error('bdry_vals must always have 6 rows');
end


% Create the right-hand side for each flux condition. Note boundary
% conditions are assumed to be the same across all tests, so they are left
% as constant.
b_allstims = inflows_list + b_bccoeff*bdry_vals;

%Solve the system (phasor at all locations) for all right-hand sides
Phi = feval(matinv_func,A,b_allstims);
%Always output from this function simulated observations
%Calculate simulated observations
sim_obs = zeros(num_obs,1);
for i = 1:1:num_obs
    stim_type = obs_list(i,1);
    obs_type = obs_list(i,2);
    sim_obs(i) = (obsweight_list(:,obs_type)).'*Phi(:,stim_type);
end

%If more than one output argument, give the full set of Phi values for all
%different stimulations (source terms). Phi is thus (num_cells x num_stims)
if nargout > 1
    varargout{1} = Phi;
end

%If more than two output arguments, the third output argument is the
%sensitivity of all simulated observations to spatially-distributed K and
%Ss throughout the model.
if nargout > 2
    
    %Sensitivity with respect to K (num_cells), then with respect to Ss
    %(num_cells), then with respect to Sy (num_x*num_y)
    H = cell(3,1);
    H{1} = zeros(num_obs,num_cells);
    H{2} = zeros(num_obs,num_cells);
    H{3} = zeros(num_obs,num_x*num_y);

    %Solve for the adjoint variable first
    lambda = feval(matinv_func,A.',obsweight_list);    
    
    %Part 1 - calculate db/ds - dA/ds*Phi where s (parameter) is K
    
    %Calculate the coefficients for all off-diagonal coefficients of dA*phi
    dAdK = zeros(num_cells,6);
    dAdK_bottom = zeros(num_y,num_x,num_z);
    for k = 2:1:(num_z)
        for j = 1:1:(num_x)
            dAdK_bottom(:,j,k) = 2*K(:,j,k-1).^2.*dy*dx(j)*dz(k)./ ...
                ((K(:,j,k)*dz(k-1) + K(:,j,k-1)*dz(k)).^2);
        end
    end
    dAdK(:,1) = reshape(dAdK_bottom,num_cells,1);
    dAdK_west = zeros(num_y,num_x,num_z);
    for j = 2:1:(num_x)
        for k = 1:1:(num_z)
            dAdK_west(:,j,k) = 2*K(:,j-1,k).^2.*dy*dx(j)*dz(k)./ ...
                ((K(:,j,k)*dx(j-1) + K(:,j-1,k)*dx(j)).^2);
        end
    end
    dAdK(:,2) = reshape(dAdK_west,num_cells,1);
    dAdK_south = zeros(num_y,num_x,num_z);
    for i = 2:1:(num_y)
        for k = 1:1:(num_z)
            dAdK_south(i,:,k) = 2*K(i-1,:,k).^2.*dx.'*dy(i)*dz(k)./ ...
                ((K(i,:,k)*dy(i-1) + K(i-1,:,k)*dy(i)).^2);
        end
    end
    dAdK(:,3) = reshape(dAdK_south,num_cells,1);
    dAdK_north = zeros(num_y,num_x,num_z);
    for i = 1:1:(num_y-1)
        for k = 1:1:(num_z)
            dAdK_north(i,:,k) = 2*K(i+1,:,k).^2.*dx.'*dy(i)*dz(k)./ ...
                ((K(i,:,k)*dy(i+1) + K(i+1,:,k)*dy(i)).^2);
        end
    end
    dAdK(:,4) = reshape(dAdK_north,num_cells,1);
    dAdK_east = zeros(num_y,num_x,num_z);
    for j = 1:1:(num_x-1)
        for k = 1:1:(num_z)
            dAdK_east(:,j,k) = 2*K(:,j+1,k).^2.*dy*dx(j)*dz(k)./ ...
                ((K(:,j,k)*dx(j+1) + K(:,j+1,k)*dx(j)).^2);
        end
    end
    dAdK(:,5) = reshape(dAdK_east,num_cells,1);
    dAdK_top = zeros(num_y,num_x,num_z);
    for k = 1:1:(num_z-1)
        for j = 1:1:(num_x)
            dAdK_top(:,j,k) = 2*K(:,j,k+1).^2.*dy*dx(j)*dz(k)./ ...
                ((K(:,j,k)*dz(k+1) + K(:,j,k+1)*dz(k)).^2);
        end
    end
    dAdK(:,6) = reshape(dAdK_top,num_cells,1);
    
    db_bccoeffs_dK = sparse(num_cells,6);
    if num_x > 1
        db_bccoeff_dK_west = zeros(num_y,num_x,num_z);
        if xtypes(1) == 1
            for k = 1:1:num_z
                db_bccoeff_dK_west(:,1,k) = 2.*dy*dz(k)/dx(1);
            end
        end
        db_bccoeffs_dK(:,1) = reshape(db_bccoeff_dK_west,num_cells,1);
        db_bccoeff_dK_east = zeros(num_y,num_x,num_z);
        if xtypes(2) == 1
            for k = 1:1:num_z
                db_bccoeff_dK_east(:,end,k) = 2.*dy*dz(k)/dx(end);
            end
        end
        db_bccoeffs_dK(:,2) = reshape(db_bccoeff_dK_east,num_cells,1);
    end
    if num_y > 1
        db_bccoeff_dK_south = zeros(num_y,num_x,num_z);
        if ytypes(1) == 1
            for k = 1:1:num_z
                db_bccoeff_dK_south(1,:,k) = 2.*dx.'*dz(k)/dy(1);
            end
        end
        db_bccoeffs_dK(:,3) = reshape(db_bccoeff_dK_south,num_cells,1);
        db_bccoeff_dK_north = zeros(num_y,num_x,num_z);
        if ytypes(2) == 1
            for k = 1:1:num_z
                db_bccoeff_dK_north(end,:,k) = 2.*dx.'*dz(k)/dy(end);
            end
        end
        db_bccoeffs_dK(:,4) = reshape(db_bccoeff_dK_north,num_cells,1);
    end
    if num_z > 1
        db_bccoeff_dK_bottom = zeros(num_y,num_x,num_z);
        if ztypes(1) == 1
            for j = 1:1:num_x
                db_bccoeff_dK_bottom(:,j,1) = 2.*dy*dx(j)/dz(1);
            end
        end
        db_bccoeffs_dK(:,5) = reshape(db_bccoeff_dK_bottom,num_cells,1);
        db_bccoeff_dK_top = zeros(num_y,num_x,num_z);
        if ztypes(2) == 1
            for j = 1:1:num_x
                db_bccoeff_dK_top(:,j,end) = 2.*dy*dx(j)/dz(end);
            end
        end
        db_bccoeffs_dK(:,6) = reshape(db_bccoeff_dK_top,num_cells,1);
    end
    dbdK_minus_dAdKPhi = cell(num_stims,1);
    for t = 1:1:num_stims
        dbdK_minus_dAdKPhi_diags = zeros(num_cells,7);
        Phi_currstim = Phi(:,t);
        bdry_vals_currstim = bdry_vals(:,t);
        Phic = reshape(Phi_currstim,num_y,num_x,num_z);
        %Create and initialize all offset-phi arrays
        Phiw = zeros(num_y,num_x,num_z);
        Phin = Phiw; Phis = Phiw; Phie = Phiw; Phit = Phiw; Phib = Phiw;
        Phis(2:end,:,:) = Phic(1:(end-1),:,:); Phis = reshape(Phis,num_cells,1);
        Phin(1:(end-1),:,:) = Phic(2:end,:,:); Phin = reshape(Phin,num_cells,1);
        Phiw(:,2:end,:) = Phic(:,1:(end-1),:); Phiw = reshape(Phiw,num_cells,1);
        Phie(:,1:(end-1),:) = Phic(:,2:end,:); Phie = reshape(Phie,num_cells,1);
        Phib(:,:,2:end) = Phic(:,:,1:(end-1)); Phib = reshape(Phib,num_cells,1);
        Phit(:,:,1:(end-1)) = Phic(:,:,2:end); Phit = reshape(Phit,num_cells,1);
        Phic = reshape(Phic,num_cells,1);
        
        %Non-main diagonals of db-dAphi (no contribution from db, so these are equal to
        %-dA*phi)
        dbdK_minus_dAdKPhi_diags(:,1) = dAdK(:,1).*(Phic - Phib);
        dbdK_minus_dAdKPhi_diags(:,2) = dAdK(:,2).*(Phic - Phiw);
        dbdK_minus_dAdKPhi_diags(:,3) = dAdK(:,3).*(Phic - Phis);
        dbdK_minus_dAdKPhi_diags(:,5) = dAdK(:,4).*(Phic - Phin);
        dbdK_minus_dAdKPhi_diags(:,6) = dAdK(:,5).*(Phic - Phie);
        dbdK_minus_dAdKPhi_diags(:,7) = dAdK(:,6).*(Phic - Phit);
        %Main diagonal of db-dAphi. db contributes to main diagonal (first
        %term), and effect of boundary conditions also affects main
        %diagonal of dA*Phi (last term)
        dbdK_minus_dAdKPhi_diags(:,4) = db_bccoeffs_dK*bdry_vals_currstim ...
            -(sum(dbdK_minus_dAdKPhi_diags,2) + db_bccoeffs_dK*ones(6,1).*Phic);
        dbdK_minus_dAdKPhi{t} = spalloc(num_cells,num_cells,7*num_cells);
        dbdK_minus_dAdKPhi{t} = spdiags(dbdK_minus_dAdKPhi_diags, ...
            [num_x*num_y num_y 1 0 -1 -num_y -num_x*num_y],...
            dbdK_minus_dAdKPhi{t});
    end

    for i = 1:1:num_obs
        test_type = obs_list(i,1);
        obs_type = obs_list(i,2);
        %Real and imaginary parts of this matrix represent the sensitivity
        %of the real and imaginary components of the phasor to K
        H{1}(i,:) = lambda(:,obs_type).'*(dbdK_minus_dAdKPhi{test_type});
    end
    
    [dxg, dyg, dzg] = meshgrid(dx,dy,dz);
    dxdydz = dxg.*dyg.*dzg;
    clear dxg dyg dzg
    
    %Part 2 - calculate db/ds - dA/ds*Phi where s (parameter) is Ss
    
    %b is not at all dependent on Ss, only sensitivity to Ss comes from
    %central diagonal of A matrix
    dbdSs_minus_dAdSsPhi = cell(num_stims,1);
    for t = 1:1:num_stims
        Phi_currstim = Phi(:,t);
        %Phic = reshape(Phi_currstim,num_y,num_x,num_z);
%         dbdSs_minus_dAdSsPhi_coldiags = omega.*1i.*reshape(dxdydz,num_cells,1)...
%             .*Phic;
        dbdSs_minus_dAdSsPhi_coldiags = -omega.*1i.*reshape(dxdydz,num_cells,1)...
            .*Phi_currstim;

        dbdSs_minus_dAdSsPhi{t} = spalloc(num_cells,num_cells,num_cells);
        dbdSs_minus_dAdSsPhi{t} = spdiags(dbdSs_minus_dAdSsPhi_coldiags,[0],...
            dbdSs_minus_dAdSsPhi{t});
    end
    
    for i = 1:1:num_obs
        test_type = obs_list(i,1);
        obs_type = obs_list(i,2);
        %Real and imaginary parts of this matrix represent the sensitivity
        %of the real and imaginary components of the phasor to Ss
        H{2}(i,:) = lambda(:,obs_type).'*dbdSs_minus_dAdSsPhi{test_type};
    end
    
    %NOTE***: In beta / not verified or supported yet!
    %Part 3 - calculate db/ds - dA/ds*Phi where s (parameter) is Sy
    dbdSy_minus_dAdSyPhi = cell(num_stims,1);
    for t = 1:1:num_stims
        Phi_currstim = Phi(:,t);
        
        %b is not at all dependent on Sy, only sensitivity to Sy comes from
        %central diagonal and below cell diagonal (num_x*num_y away).
        dAdSy_center = zeros(num_y,num_x,num_z);
        dAdSy_bottom = zeros(num_y,num_x,num_z);
        %TODO - Verify. Make sure it matches with by-hand
        %computation of inner product with matrix
        if ztypes(2) == 2
            for i = 1:1:num_y
                for j = 1:1:num_x
                    dAdSy_center(i,j,end) = omega*1i*dx(j)*dy(i)*(1+(dz(end)/(dz(end) + dz(end-1))));
                    dAdSy_bottom(i,j,end) = -omega*1i*dx(j)*dy(i)*(dz(end)/(dz(end) + dz(end-1)));
                end
            end
        end
        Phic = reshape(Phi_currstim,num_y,num_x,num_z);
        %Phib = zeros(num_y,num_x,num_z);
        %Phib(:,:,2:end) = Phic(:,:,1:(end-1)); Phib = reshape(Phib,num_cells,1);
        dbdSy_minus_dAdSyPhi{t} = spalloc(num_cells,num_cells,num_cells);
        dbdSy_minus_dAdSyPhi_coldiags = reshape(dAdSy_center.*Phic + dAdSy_bottom.*Phic,num_cells,1);
        dbdSy_minus_dAdSyPhi{t} = spdiags(dbdSy_minus_dAdSyPhi_coldiags,[0],...
            dbdSy_minus_dAdSyPhi{t});
    end
    
    for i = 1:1:num_obs
        test_type = obs_list(i,1);
        obs_type = obs_list(i,2);
        %Real and imaginary parts of this matrix represent the sensitivity
        %of the real and imaginary components of the phasor to Sy
        %TODO - Verify. Hack for integrating over only the top
        %layer.
        Sy_sens_mat = reshape(full(lambda(:,obs_type).'*(dbdSy_minus_dAdSyPhi{test_type})),num_y,num_x,num_z);
        H{3}(i,:) = reshape(Sy_sens_mat(:,:,end),num_y*num_x,1);
    end
    
    varargout{2} = H;
end