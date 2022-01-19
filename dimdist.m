function distmat = dimdist(X,varargin)

%dimdist: Function for creating a set of matrices giving the distance along
%each dimension
%
%[distmat] = dimdist(X,{Y})
%
%where:
%   -distmat is a (dim x 1) cell array where each cell contains a matrix of
%   distances between the points along one dimension. The dimension of each
%   matrix will be numpoints1 x numpoints2
%   -X is a (numpoints1 x dim) matrix where each row contains the location
%   of a point
%   -Y (optional) is a (numpoints2 x dim) matrix containing a second set of
%   points. If Y is not supplied, Y is set to be equal to X.
%
% © 2012-2013 Michael Cardiff, Warren Barrash, and Peter K. Kitanidis, All
% Rights Reserved.


if nargin == 2 && ~isempty(varargin{1})
    Y = varargin{1};
else
    Y = X;
end

if nargin == 3
    signed = varargin{2};
else
    signed = 0;
end

m = size(X,1);
n = size(Y,1);
dim = size(X,2);

if size(Y,2) ~= dim
    error('Both vectors must have the same number of dimensions for coordinate definitions');
end

distmat = zeros(m,n,dim);
for i = 1:1:dim
    v1 = X(:,i); v2 = Y(:,i);
    [m1,m2] = meshgrid(v2,v1);
    if signed == 0
        distmat(:,:,i) = abs(m1-m2);
    else
        distmat(:,:,i) = m1-m2;
    end
end