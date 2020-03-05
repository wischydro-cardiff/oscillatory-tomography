function weights = grid_idw(interp_locs,grid_x,varargin)

% grid_idw: Function that produces inverse distance weighting coefficients
% given a set of grid nodes and locations at which values are to be
% interpolated. 
%
% Note: This function is built only for interpolation where the basis is a
% rectangular grid and the neighborhood of interpolation is, at most, the
% set of maximum 2^dim closest points. More general functions must be used
% for non-rectangular grids.
%
% Syntax:
% [weights] = grid_idw(interp_locs,grid_x,[grid_y],[grid_z])
% where:
%   -weights is a (num_interpolants x num_gridpts) matrix of the
%   interpolation weights. num_interpolants is the number of points for
%   which interpolation is being performed, and num_gridpts is the number
%   of points in the grid. NOTE: The ordering of the columns values in
%   weights follows the convention of meshgrid, i.e. it iterates through y
%   first, then through x, then through z (if 3D).
%   -interp_locs is a (num_interpolants x dim) matrix of locations for
%   which interpolation is to be carried out, where dim is the dimension of
%   the problem
%   -grid_x is a (num_x x 1) monotonically increasing vector of x
%   coordinates where num_x is the number of unique x locations in the grid
%   -grid_y (optional) is a (num_y x 1) monotonically increasing vector of
%   y coordinates where num_y is the number of unique y locations in the
%   grid
%   -grid_z (optional) is a (num_z x 1) monotonically increasing vector of
%   z coordinates where num_z is the number of unique z locations in the
%   grid
%
% TODO: Add documentation / change documentation Added option to roughly
% approximate partially-penetrating conditions if a beginning and end z
% location are given. Weights are equally distributed between all layers
% that have some contribution between the z maximum and z minimum.

grid_y = [];
grid_z = [];

[grid_y grid_z] = process_extra_args(varargin,grid_y,grid_z);

%Error checking on coordinates supplied
if isempty(grid_y) && ~isempty(grid_z)
    error('Coordinates must be given in order. No intermediate grid coordinate list may be empty.');
end

%Figure out the dimension of the problem being solved
if isempty(grid_y)
    dim = 1;
    numx = numel(grid_x);
elseif isempty(grid_z)
    dim = 2;
    numx = numel(grid_x);
    numy = numel(grid_y);
else
    dim = 3;
    numx = numel(grid_x);
    numy = numel(grid_y);
    numz = numel(grid_z);
end

ln_src = 0;
if size(interp_locs,2) == 4
    ln_src = 1;
end

if size(interp_locs,2) < dim
    error('Dimensions of interp_locs appears to be too small');
end

if (size(interp_locs,2) > dim) && (dim ~= 3)
    error('Size of interp_locs appears to be too large.');
end

%Check that interpolants are within the grid, and check that all grid
%coordinate vectors are monotonically increasing
if (max(interp_locs(:,1)) > max(grid_x)) || (min(interp_locs(:,1)) < min(grid_x))
    error('X locations to be interpolated do not fall within the span of the grid');
end
if ~all((grid_x(2:end) - grid_x(1:(end-1))) > 0)
    error('Grid X locations must be a monotonically increasing vector of values');
end
if dim >= 2
    if (max(interp_locs(:,2)) > max(grid_y)) || (min(interp_locs(:,2)) < min(grid_y))
        error('Y locations to be interpolated do not fall within the span of the grid');
    end
    if ~all((grid_y(2:end) - grid_y(1:(end-1))) > 0)
        error('Grid Y locations must be a monotonically increasing vector of values');
    end
end
if dim >= 3
    if (max(interp_locs(:,3)) > max(grid_z)) || (min(interp_locs(:,3)) < min(grid_z))
        error('Z locations to be interpolated do not fall within the span of the grid');
    end
    if ln_src == 1
        if (max(interp_locs(:,4)) > max(grid_z)) || (min(interp_locs(:,4)) < min(grid_z))
            error('Z locations to be interpolated do not fall within the span of the grid');
        end
    end
    if ~all((grid_z(2:end) - grid_z(1:(end-1))) > 0)
        error('Grid Z locations must be a monotonically increasing vector of values');
    end
end

num_interpolants = size(interp_locs,1);
num_nearest_max = 2^dim;
if dim == 1
    num_gridpts = numx;
elseif dim == 2
    num_gridpts = numx*numy;
elseif dim == 3
    num_gridpts = numx*numy*numz;
end

weights = spalloc(num_interpolants,num_gridpts,num_interpolants*num_nearest_max);
if dim == 1
    for i = 1:1:num_interpolants
        %Pre-allocated neighbor_list and invdist_list to their maximum size
        neighbor_list = zeros(num_nearest_max,dim);
        invdist_list = zeros(num_nearest_max,1);
        
        %Find the grid cell just beyond the interpolation point
        upx_ind = find(grid_x >= interp_locs(i,1),1);
        %If the interp point falls right on a value in a particular
        %dimension, only use that point. Otherwise, use the current point
        %(which is above) and another point below.
        %Note: This also protects against the case where the interpolation
        %point is exactly at a boundary.
        if (grid_x(upx_ind) == interp_locs(i,1))
            xn_vec = [upx_ind];
        else
            xn_vec = [upx_ind-1, upx_ind];
        end
        %Create and count list of all neighbors
        num_neighbors = numel(xn_vec);
        neighbor_list([1:num_neighbors],:) = (combvec(xn_vec))';
        neighbor_list = neighbor_list([1:num_neighbors],:);
        invdist_list = zeros(num_neighbors,1);
        %Calculate inverse distances and weights for all neighboring points
        invdist_list = (grid_x(neighbor_list(:,1)) - interp_locs(i,1)).^2;
        invdist_list = 1./(invdist_list.^.5);
        invdist_list(invdist_list == inf) = 1./eps;
        weight_list = invdist_list./(sum(invdist_list));
        grid_inds = neighbor_list(:,1);
        weights(i,grid_inds) = weight_list;
    end
elseif dim == 2
    for i = 1:1:num_interpolants
        %Pre-allocated neighbor_list and invdist_list to their maximum size
        neighbor_list = zeros(num_nearest_max,dim);
        invdist_list = zeros(num_nearest_max,1);
        
        %Find the grid cell just beyond the interpolation point
        upx_ind = find(grid_x >= interp_locs(i,1),1);
        upy_ind = find(grid_y >= interp_locs(i,2),1);
        %If the interp point falls right on a value in a particular
        %dimension, only use that point. Otherwise, use the current point
        %(which is above) and another point below.
        %Note: This also protects against the case where the interpolation
        %point is exactly at a boundary.
        if (grid_x(upx_ind) == interp_locs(i,1))
            xn_vec = [upx_ind];
        else
            xn_vec = [upx_ind-1, upx_ind];
        end
        if (grid_y(upy_ind) == interp_locs(i,2))
            yn_vec = [upy_ind];
        else
            yn_vec = [upy_ind-1, upy_ind];
        end
        %Create and count list of all neighbors
        num_neighbors = numel(xn_vec)*numel(yn_vec);
        neighbor_list([1:num_neighbors],:) = (combvec(xn_vec,yn_vec))';
        neighbor_list = neighbor_list([1:num_neighbors],:);
        invdist_list = zeros(num_neighbors,1);
        %Calculate inverse distances and weights for all neighboring points
        invdist_list = (grid_x(neighbor_list(:,1)) - interp_locs(i,1)).^2 + ...
            (grid_y(neighbor_list(:,2)) - interp_locs(i,2)).^2;
        invdist_list = 1./(invdist_list.^.5);
        invdist_list(invdist_list == inf) = 1./eps;
        weight_list = invdist_list./(sum(invdist_list));
        grid_inds = (neighbor_list(:,1)-1)*numy + neighbor_list(:,2);
        weights(i,grid_inds) = weight_list;
    end
elseif dim == 3
    if ln_src == 0
        for i = 1:1:num_interpolants
            %Pre-allocated neighbor_list and invdist_list to their maximum size
            neighbor_list = zeros(num_nearest_max,dim);
            invdist_list = zeros(num_nearest_max,1);
            
            %Find the grid cell just beyond the interpolation point
            upx_ind = find(grid_x >= interp_locs(i,1),1);
            upy_ind = find(grid_y >= interp_locs(i,2),1);
            upz_ind = find(grid_z >= interp_locs(i,3),1);
            %If the interp point falls right on a value in a particular
            %dimension, only use that point. Otherwise, use the current point
            %(which is above) and another point below.
            %Note: This also protects against the case where the interpolation
            %point is exactly at a boundary.
            if (grid_x(upx_ind) == interp_locs(i,1))
                xn_vec = [upx_ind];
            else
                xn_vec = [upx_ind-1, upx_ind];
            end
            if (grid_y(upy_ind) == interp_locs(i,2))
                yn_vec = [upy_ind];
            else
                yn_vec = [upy_ind-1, upy_ind];
            end
            if (grid_z(upz_ind) == interp_locs(i,3))
                zn_vec = [upz_ind];
            else
                zn_vec = [upz_ind-1, upz_ind];
            end
            %Create and count list of all neighbors
            num_neighbors = numel(xn_vec)*numel(yn_vec)*numel(zn_vec);
            neighbor_list([1:num_neighbors],:) = (combvec(xn_vec,yn_vec,zn_vec))';
            neighbor_list = neighbor_list([1:num_neighbors],:);
            invdist_list = zeros(num_neighbors,1);
            %Calculate inverse distances and weights for all neighboring points
            invdist_list = (grid_x(neighbor_list(:,1)) - interp_locs(i,1)).^2 + ...
                (grid_y(neighbor_list(:,2)) - interp_locs(i,2)).^2 + ...
                (grid_z(neighbor_list(:,3)) - interp_locs(i,3)).^2;
            invdist_list = 1./(invdist_list.^.5);
            invdist_list(invdist_list == inf) = 1./eps;
            weight_list = invdist_list./(sum(invdist_list));
            grid_inds = (neighbor_list(:,3)-1)*numy*numx + ...
                (neighbor_list(:,1)-1)*numy + neighbor_list(:,2);
            weights(i,grid_inds) = weight_list;
        end
    else
        for i = 1:1:num_interpolants
            %Pre-allocated neighbor_list and invdist_list to their maximum size
            neighbor_list = zeros(num_nearest_max,dim);
            invdist_list = zeros(num_nearest_max,1);
            
            %Find the grid cell just beyond the interpolation point
            upx_ind = find(grid_x >= interp_locs(i,1),1);
            upy_ind = find(grid_y >= interp_locs(i,2),1);
            upz_ind = find(grid_z >= interp_locs(i,4),1);
            dnz_ind = find(grid_z <= interp_locs(i,3),1,'last');
            %If the interp point falls right on a value in a particular
            %dimension, only use that point. Otherwise, use the current point
            %(which is above) and another point below.
            %Note: This also protects against the case where the interpolation
            %point is exactly at a boundary.
            if (grid_x(upx_ind) == interp_locs(i,1))
                xn_vec = [upx_ind];
            else
                xn_vec = [upx_ind-1, upx_ind];
            end
            if (grid_y(upy_ind) == interp_locs(i,2))
                yn_vec = [upy_ind];
            else
                yn_vec = [upy_ind-1, upy_ind];
            end
            zn_vec = dnz_ind:upz_ind;
            %Create and count list of all neighbors
            num_neighbors = numel(xn_vec)*numel(yn_vec)*numel(zn_vec);
            neighbor_list([1:num_neighbors],:) = (combvec(xn_vec,yn_vec,zn_vec))';
            neighbor_list = neighbor_list([1:num_neighbors],:);
            invdist_list = zeros(num_neighbors,1);
            %Calculate inverse distances and weights for all neighboring points
            invdist_list = (grid_x(neighbor_list(:,1)) - interp_locs(i,1)).^2 + ...
                (grid_y(neighbor_list(:,2)) - interp_locs(i,2)).^2;
            invdist_list = 1./(invdist_list.^.5);
            invdist_list(invdist_list == inf) = 1./eps;
            weight_list = invdist_list./(sum(invdist_list));
            %Modification for length of line contained per cell
            %TODO -Insert here
            grid_inds = (neighbor_list(:,3)-1)*numy*numx + ...
                (neighbor_list(:,1)-1)*numy + neighbor_list(:,2);
            weights(i,grid_inds) = weight_list;
        end
    end
end