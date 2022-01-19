function [coords, varargout] = plaid_cellcenter_coord(domain)

%plaid_coord: Function to create a list of the center coordinates of points
%in a plaid grid as a vector of 2- or 3-dimensional coordinates.
%
%[coords, {grid_pts}] = plaid_coord(domain)
%
%where:
%   OUTPUTS:
%   -coords is a numcells x 3 matrix where each row contains the
%   coordinates of one point in the plaid grid. numcells = numx * numy (*
%   numz, if supplied)
%   -grid_pts is a 3 x 1 cell array, where each cell contains a matrix
%   (numy x numx x numz) containing the spatial coordinates for one
%   dimension in the plaid grid. grid_pts is exactly the same output as
%   produced by meshgrid.
%   INPUTS:
%   -domain is a structure with fields x,y, and z each containing the cell
%   boundaries along this dimension
%
% Code by Michael Cardiff
% 10/2014, last updated 1/2016

if ~isfield(domain,'x')
    error('X boundary information must be supplied in domain.x');
else
    xb = domain.x;
end

if ~isfield(domain,'y')
    error('Y boundary information must be supplied in domain.y');
else
    yb = domain.y;
end

if ~isfield(domain,'z')
    error('Z boundary information must be supplied in domain.z');
else
    zb = domain.z;
end


xc = (xb(2:end) + xb(1:(end-1)))./2;
yc = (yb(2:end) + yb(1:(end-1)))./2;
zc = (zb(2:end) + zb(1:(end-1)))./2;

num_x = numel(xc);
num_y = numel(yc);
num_z = numel(zc);

if num_z > 1
    num_pts = num_x*num_y*num_z;
    [xg, yg, zg] = meshgrid(xc,yc,zc);
    coords = [reshape(xg,num_pts,1) reshape(yg,num_pts,1) reshape(zg,num_pts,1)];    
    if nargout == 2
        varargout{1} = {xg yg zg};
    end
else
    [xg, yg] = meshgrid(xc,yc);
    num_pts = num_x*num_y;
    coords = [reshape(xg,num_pts,1) reshape(yg,num_pts,1)];
    if nargout == 2
        varargout{1} = {xg yg};
    end
end


