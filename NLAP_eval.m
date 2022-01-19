function [NLAP] = NLAP_eval(y,X,s,beta,Q,R,h_function)

%negloglike_eval:Function which calculates the negative log-likelihood
%using a general model (h_function), and given either a Q matrix or a
%function that calculates a Q matrix-vector product
%
% [negloglike] = negloglike_eval(y,X,s,beta,Q,R,h_function)
% where:
%   OUTPUTS:
%       -NLAP is the negative log a-posteriori objective function value
%       (consisting of the data misfit contribution and the "prior" contribution)
%   INPUTS:
%       -y is the (m x 1) vector of data
%       -X is the (n x p) matrix with each column containing the values of
%       the drift functions
%       -s is the estimate of the parameters (n x 1)
%       -beta is the estimate of the drift coefficients (p x 1);
%       -Q is the parameter covariance matrix (n x n), or an anonymous
%       function that takes a vector v as input and calculates the product
%       Q*v.
%       -R is the data error covariance matrix (m x m)
%       -h_function is an anonymous function that takes the parameter
%       estimates vector as inputs and calculates the simulated data vector
%       (n --> m)
%
% Note: If Q is an anonymous function, Q^-1 is approximated in a matrix
% multiplication step using minres with a specified tolerance and number of
% iterations. In the future, these tolerances and iteration parameters will
% be added to the list of possible inputs for flexibility.
%
% Code by Michael Cardiff, 2010-2016

resid = y - h_function(s);

NLAP_resid_part = 0.5.*resid'*(R\resid);
xi = s-X*beta;

if isa(Q,'function_handle')
    NLAP_geostat_part = 0.5.*xi'*minres(Q,xi,1e-4,10000);
else
    NLAP_geostat_part = 0.5.*xi'*(Q\xi);
end
    
NLAP = NLAP_resid_part + NLAP_geostat_part;
