function distmat = eucdist(X,varargin)

%eucdist: Modified function for calculating distance between a set of
%points and organizing data in a distance matrix
%
% © 2012-2013 Michael Cardiff, Warren Barrash, and Peter K. Kitanidis, All
% Rights Reserved.


if nargin == 2
    Y = varargin{1};
else
    Y = X;
end

m = size(X,1);
n = size(Y,1);
dim = size(X,2);

if size(Y,2) ~= dim
    error('Both vectors must have the same number of dimensions for coordinate definitions');
end

distmat = zeros(m,n);
for i = 1:1:dim
    v1 = X(:,i); v2 = Y(:,i);
    [m1,m2] = meshgrid(v2,v1);
    distmat = distmat + (m1-m2).^2;
end
distmat = distmat.^.5;