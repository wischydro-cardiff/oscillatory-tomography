function [delta_resid, standard_error, orthonorm_resid, cR, Q2] = geostat_resid_compute(y,H,X,R,Q)

%geostat_resid_compute: Function which computes orthonormal residuals and
%other statistics for a geostatistically-based inverse problem. Assumes
%LINEAR system, but can be called for linearized quasi-linear systems.
%   [delta_resid, standard_error, orthonorm_resid, cR, Q2] = 
%     geostat_resid_compute(y,H,X,R,Q)
%   where:
%       -y is the data (linearized, as appropriate)
%       -H is the forward model, or sensitivity matrix in the linearized
%       case
%       -X is the assignment / drift matrix
%       -R is the error covariance matrix
%       -Q is the spatial covariance matrix of the parameters

m = size(y,1);
n = size(Q,1);
p = size(X,2);

PSI = H*Q*H'+R;
PHI = H*X;

P = null(PHI')';
Pyy = P'*inv(P*PSI*P')*P;
T = (orth(Pyy))';
delta_resid = T*y;
vars = diag(T*PSI*T');
standard_error = sqrt(vars);
orthonorm_resid = delta_resid./standard_error;

Q2 = 1/(m-p).*sum(orthonorm_resid.^2);
cR = Q2*exp(1/(m-p).*sum(log(vars)));
