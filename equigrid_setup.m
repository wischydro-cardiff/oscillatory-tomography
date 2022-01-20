function [domain] = equigrid_setup(x_disc, y_disc, varargin)

%equigrid_setup: This function takes simple user inputs for an equispaced
%grid and generates the "domain" structure needed, consisting of cell boundary locations.
%
%[domain] = equigrid_setup(x_disc,y_disc,[z_disc])
%
%where:
%   OUTPUTS:
%   -domain is a structure containing fields x,y, and z which each contain
%   a vector of cell boundary locations
%   INPUTS:
%   -x_disc is a 3 x 1 vector containing information about the
%   x-discretization. The three elements are the minimum x value, maximum x
%   value, and the number of cells in the x dimension.
%   -y_disc is a 3 x 1 vector for the y discretization, similar to x_disc
%   -z_disc (optional) is a 3 x 1 vector for the z discretization
%
% Code by Michael Cardiff
% 10/2014, Last Updated: 1/2016 

num_baseargs = 2;
num_inargs = nargin;

%Setup discretization parameters, assuming equal spacing along a dimension.
xbs = [x_disc(1) x_disc(2)];
num_x = x_disc(3);
ybs = [y_disc(1) y_disc(2)];
num_y = y_disc(3);

%If Z discretization is given use it... otherwise assume unit thickness and
%one cell.
if num_inargs > num_baseargs
    z_disc = varargin{1};
else
    z_disc = [0 1 1];
end
zbs = [z_disc(1) z_disc(2)];
num_z = z_disc(3);

%Cell boundaries - vectors
xb = linspace(xbs(1),xbs(2),num_x+1);
yb = linspace(ybs(1),ybs(2),num_y+1);
zb = linspace(zbs(1),zbs(2),num_z+1);

domain = struct('x',xb,'y',yb,'z',zb);
