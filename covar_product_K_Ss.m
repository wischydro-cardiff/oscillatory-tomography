function [Qproduct] = covar_product_K_Ss(QK_row,QSs_row,invec,numx,numy,varargin)

%covar_product_K_Ss: Function for computing product of a vector with the
%covariance matrix for the aquifer parameters K and Ss. Code assumes that
%the covariance matrix Q has the structure Q = [QK zeros(n,n); zeros(n,n)
%QSs], where QK and QSs are determined by a spatial covariance function,
%i.e. the variogram. This means no correlation is assumed between K values
%and Ss values within the aquifer. This code utilizes the toeplitz approach
%to matrix-vector products, thus only the first rows of the Q matrices need
%to be supplied
%
%[Qproduct] = covar_product_K_Ss(QK_row,QSs_row,invec,numx,numy,{numz})
%
%Where:
%   -Qproduct is the product of the Q matrix with invec 
%   -QK_row (numcells x 1) is the first row of the covariance matrix for K
%   values, where numcells is the number of cells in the model.
%   -QSs_row (numcells x 1) is the first row of the covariance matrix for
%   Ss values
%   -invec ((numcells*2) x 1) is the vector being multiplied by Q
%   -numx is the number of cells in the x direction
%   -numy is the number of cells in the y direction
%   -numz (optional) is the number of cells in the z direction
%
% Code by Michael Cardiff
% 10/2014, Last Updated: 12/2014 


num_reqin = 5;
dim = 2;

if nargin > num_reqin
    numz = varargin{1};
    dim = 3;
end

numK = numel(QK_row);
numSs = numel(QSs_row);
numtot = numK+numSs;

%"invec" must contain values corresponding to all parameter values
%for K followed by all parameter values for Ss, in order.
invec_Kpart = invec(1:numK);
invec_Sspart = invec((numK+1):(numK+numSs));

if dim == 2
    QK_prod = toepmat_vector_math(QK_row,'*',invec_Kpart,2,[numy numx]);
    QSs_prod = toepmat_vector_math(QSs_row,'*',invec_Sspart,2,[numy numx]);
else
    QK_prod = toepmat_vector_math(QK_row,'*',invec_Kpart,3,[numy numx numz]);
    QSs_prod = toepmat_vector_math(QSs_row,'*',invec_Sspart,3,[numy numx numz]);
end

Qproduct = [QK_prod; QSs_prod];
   