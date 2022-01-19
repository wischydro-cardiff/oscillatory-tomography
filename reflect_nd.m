function [reflected_point] = reflect_nd(point_reflectme,plane_vec)

%reflect_nd: Function which reflects a point over a plane, in n dimensions.
%Works by first translating problem so that plane passes through origin,
%then performs reflection using simple rule, then back translates.
%syntax:
%[reflected_point] = reflect_nd(point_reflectme,plane_vec)
%where:
%   -reflected_point is a (d x 1) vector containing the coordinates of the
%   reflected point, where d is the dimension of the problem
%   -point_reflectme is a (d x 1) vector giving the point's coordinates
%   -plane_vec is a (d+1 x 1) vector which defines the plane. The product
%   plane_vec'*([X; 1]) = 0 should give the formula for the plane. For
%   example, in 2D if plane_vec = [a; b; c], this implies ax + by + c = 0

dim = numel(plane_vec) - 1;
dim2 = numel(point_reflectme);

if size(plane_vec,2) > 1
    plane_vec = plane_vec';
end
if size(point_reflectme,2) > 1
    point_reflectme = point_reflectme';
end

plane_normal = plane_vec(1:end-1);

if dim ~= dim2
    error('Dimensions do not match');
end

if sum(abs(plane_normal)) == 0
    error('plane_vec must have at least one non-zero coefficient');
end

coord_toedit = find(plane_vec,1);

translate_amt = -plane_vec(end)/plane_vec(coord_toedit);

point_reflectme(coord_toedit) = point_reflectme(coord_toedit) - translate_amt;

%Formula from Wikipedia for reflection
%http://en.wikipedia.org/wiki/Reflection_%28mathematics%29
reflected_point = point_reflectme - 2*(point_reflectme'*plane_normal)/(plane_normal'*plane_normal)*plane_normal;

reflected_point(coord_toedit) = reflected_point(coord_toedit) + translate_amt;

