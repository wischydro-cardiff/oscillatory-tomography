function [coords, varargout] = plaid_coord(x,y,varargin);

%plaid_coord: Function to create a list of the coordinates of points
%in a plaid grid as a vector of 2- or 3-dimensional coordinates.
%
%[coords, {grid_pts}] = plaid_coord(x,y,{z})
%
%where:
%   -coords is a numcells x 3 matrix where each row contains the
%   coordinates of one point in the plaid grid. numcells = numx * numy (*
%   numz, if supplied)
%   -grid_pts is a 3 x 1 cell array, where each cell contains a matrix
%   (numy x numx x numz) containing the spatial coordinates for one
%   dimension in the plaid grid. grid_pts is exactly the same output as
%   produced by meshgrid.
%   -x (numx x 1) is the vector of points in the x dimension
%   -y (numy x 1) is the vector of points in the y dimension
%   -z (numz x 1, optional) is the vector of points in the z dimension
%
% Code by Michael Cardiff
% 10/2014, last updated 12/2014

xct = numel(x);
yct = numel(y);

if nargin == 3
    z = varargin{1};
    zct = numel(z);
    num_pts = xct*yct*zct;
    [xg, yg, zg] = meshgrid(x,y,z);
    coords = [reshape(xg,num_pts,1) reshape(yg,num_pts,1) reshape(zg,num_pts,1)];    
    if nargout == 2
        varargout{1} = {xg yg zg};
    end
elseif nargin == 2
    [xg, yg] = meshgrid(x,y);
    num_pts = xct*yct;
    coords = [reshape(xg,num_pts,1) reshape(yg,num_pts,1)];
    if nargout == 2
        varargout{1} = {xg yg};
    end
else
    error('Function plaid_coord only accepts 2 or 3 arguments');
end


