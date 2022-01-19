function [Q] = compute_Q_nD(unk_coords, GCF_func, GCF_params)

%compute_Q_nD.m: Function used to compute the generalized covariance matrix
%based on the given generalized covariance function and the locations of
%the unknowns. This code works for all dimensions, from 1 to N. Each row of
%the input coordinates matrix should contain the coordinates of one
%unknown. This code only deals with isotropic covariance functions.
%   [Q] = compute_Q_nD(unk_coords, GCF_func, GCF_params) 
%   where:
%       -unk_coords is a matrix containing the coordinates of each point
%       (one per row).
%       -GCF_func is a function that accepts two inputs (distance,
%       GCF parameters) and computes the GCF value for those two
%       points. Function should be vectorized so that it can handle a
%       general distance matrix.
%       -GCF_params are the GCF parameters for the given generalized
%       covariance function.

dims = size(unk_coords,2);
x = unk_coords;
dist = squareform(pdist(x));
Q = GCF_func(dist, GCF_params);