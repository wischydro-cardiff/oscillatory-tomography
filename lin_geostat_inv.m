function [s_hat, xi, beta,varargout] = lin_geostat_inv(y,X,R,Q,H)

%lin_geostat_inv.m: Function used to compute the best estimate of
%parameters given data for a linear forward model.
%   [s_hat, xi, beta,[s_postcov]] = lin_invert_solve(y,X,R,Q,H)
%   where:
%   OUTPUTS:
%       -s_hat is the best estimate of the parameters (n x 1)
%       -xi is the set of coefficients for the "basis functions" as
%       described in Snodgrass and Kitanidis (multipliers for the matrix
%       Q*H')
%       -beta is the set of drift function coefficients (p x 1)
%       -s_postcov (optional) is an n x n matrix containing the estimated
%       posterior covariance matrix
%   INPUTS:
%       -y is the data (m x 1). If solving a nonlinear case, y should be
%       passed as y-h(s_tilde) + H*s_tilde
%       -X is the assignment or drift matrix (n x p)
%       -R is the data error covariance matrix (m x m)
%       -Q is the spatial covariance matrix of the parameters (n x n), OR
%       an anonymous function that evaluates the product Q*v when a vector
%       v is supplied as the input.
%       -H is the forward model (or sensitivity matrix, in the nonlinear
%       case), an (m x n) matrix
%
%
% Note: At this point, the posterior covariance matrix will only be
% returned correctly if Q is supplied as an actual matrix. Support for
% posterior covariances when Q is supplied as an anonymous function is
% coming shortly.
%
% Note: Currently, all matrix inverse operations are performed using the
% MATLAB \ operator. Future versions will include options for using other
% approximate matrix inversion approaches.
%
% Code by Michael Cardiff, 2009-2016

num_reqout = 3;

m = size(y,1);
n = size(X,1);
p = size(X,2);

if isa(Q,'function_handle')
    %In this case, QHt should be calculated first through matrix
    %multiplication (reduces size to n x m), and then multiplied by H
    %(reduces size to m x m).
    QHt = zeros(n,m);
    for i = 1:1:m
        QHt(:,i) = Q(H(i,:)');
    end
    PSI = [H*QHt + R];  
    
elseif ismatrix(Q)
    %In this case, Q is small enough to be stored. All of the matrix-vector
    %products here should be smaller, for under-determined inverse
    %problems. Thus little need to worry about memory.
    PSI = H*Q*H' + R;
        
end

PHI = H*X;

xi_beta_vector = ([PSI, PHI; PHI', zeros(p,p)])\...
    [y; zeros(p,1)];

A = ([PSI, PHI; PHI', zeros(p,p)]);

xi = xi_beta_vector(1:m);
beta = xi_beta_vector(m+1:m+p);

if isa(Q,'function_handle')
    s_hat = X*beta + QHt*xi;
else
    s_hat = X*beta + Q*H'*xi;
end

if nargout > num_reqout
    s_postcov = Q - [H*Q; X']'*(A\[H*Q; X']);
    varargout{1} = s_postcov;
end

end