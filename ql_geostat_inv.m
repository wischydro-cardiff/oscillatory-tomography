function [s_best, beta_best, H_local, NLAP] = ql_geostat_inv(y,s_init,beta_init,X,R,Q,h_function,H_tilde_func, varargin)

%ql_geostat_inv.m: Inverse problem solver using Kitanidis' Quasi-linear
%geostatistical method. 
% Syntax:
% [s_best, beta_best, H_local, NLAP] =
% ql_geostat_inv(y,s_init,X,R,Q,forward_function,sensitivity_function,...)
% Note: In everything below:
% m = number of observations
% n = number of parameter values
% p = number of parameters in geostatistical trend
%
% where:
%   OUTPUTS:
%   -s_best is an (n x 1) list of best estimates for the parameter values
%   -beta_best is a (p x 1) matrix containing the optimal drift coefficients
%   -H_local is a local (linearized) jacobian (m x n), as evaluated at the
%   estimate s_best
%   -NLAP is the negative log-a posteriori value (objective function value)
%   obtained at the parameters s_best.
%   INPUTS
%   -s_init is an (n x 1) list of initial parameter value guesses
%   -y is an (m x 1) list of observed values
%   -Q is an (n x n) prior covariance matrix for the parameter values
%   -R is an (m x m) matrix of covariance of the measurement errors
%   -X is an (n x p) matrix, with each column representing the values of a
%   drift function.
%   -forward_function is an anonymous function with one vector input argument
%   (input of parameter values), that returns a vector (value of all
%   observations)
%   -sensitivity_function is an anonymous function with one vector input
%   argument (input of all parameter values), that returns a Jacobian
%   (sensitivity of all observations to all parameters), an (m x n) matrix
%
% Any other variables (i.e., those necessary for evaluation of
% sensitivity_function or forward_function) can be passed after the
% sensitivity function argument as follows:
% ...,'fem',fem,'max_linsearch',10,...
%
% Arguments used by ql_geostat_inv which can be passed in this section:
%    -max_linsearch: Maximum number of steps in linesearch optimization
%    -max_gradevals: Maximum number of gradient evaluations (maximum number
%    of iterations, in essence).
%    -tol_objfunc: Relative tolerance in the objective function change
%    (percent decrease)
%    -tol_s: Relative tolerance in the change in the value of s (root mean
%    square difference from last iteration).
%    -tol_linsearch: Relative tolerance in the change of the objective
%    function during the linesearch optimization (if used)
%
% NOTE: For debugging purposes, this function produces two global
% variables, s_hat and H_tilde that represent the current values of the
% parameter estimates and the current Jacobian / sensitivity matrix. This
% can be used to restart inversions, for example, if the global variables
% are addressed externally.
%
% Code by Michael Cardiff, 2010-2016

%DEBUGGING: For debugging / crashes
%Create global variables - in case function crashes due to memory or other
%errors, these variables can still be accessed outside of this program.
% global H_tilde
% global s_hat

start_time = clock;

%Default values for all optimization constants
max_linsearch = 20;
max_gradevals = 30;
tol_objfunc = .001;
tol_s = .001;
tol_linsearch = .001;
m = size(y,1);
n = size(s_init,1);
p = size(X,2);
s_tilde = s_init;
s_hat = s_tilde;
beta_tilde = beta_init;
beta_hat = beta_tilde;
nreqin = 8;

%Load auxilliary variables
num_additional = nargin - nreqin;
if num_additional > 0
    if mod(num_additional,2) == 0
        for i = 1:2:num_additional
            eval([varargin{i}, ' = varargin{', num2str(i+1), '};']);
        end
    else
        error('Wrong number of additional arguments');
    end
end

NLAP_func = @(s,beta) NLAP_eval(y,X,s,beta,Q,R,h_function);

NLAP = NLAP_func(s_tilde,beta_tilde);
NLAP_new = NLAP;

objfunc_pct = (NLAP - NLAP_new)/NLAP
s_pct = (sum(abs(s_tilde - s_hat)./s_tilde))/size(s_tilde,1)

num_gradevals = 0;
while ((objfunc_pct > tol_objfunc) && (s_pct > tol_s) && (num_gradevals < max_gradevals)) || (num_gradevals == 0)
    
    format short
    elapsed = etime(clock, start_time);

    s_tilde = s_hat;

    NLAP = NLAP_new;
    disp('Iteration         Time Elapsed        NLAP');
    disp([num2str(num_gradevals), '      ', num2str(elapsed), '       ', num2str(NLAP)]);

    H_tilde = H_tilde_func(s_tilde);
    num_gradevals = num_gradevals + 1;

    h_of_s_tilde = h_function(s_tilde);
    [s_hat, ~, beta_hat] = lin_geostat_inv((y-h_of_s_tilde+H_tilde*s_tilde),X,R,Q,H_tilde);

    NLAP_new = NLAP_func(s_hat,beta_hat);

    if max_linsearch > 0

        clear linsearch_func
        linsearch_func = @(x) NLAP_func(s_tilde + (s_hat - s_tilde).*x,...
            beta_tilde + (beta_hat - beta_tilde).*x);
        if NLAP_new < NLAP
            disp('Performing linesearch, starting at s_hat');
            lin_options = optimset('MaxIter', max_linsearch,'Display','iter','TolFun',tol_linsearch*NLAP_new);
            step = fminsearch(linsearch_func,1,lin_options)
        else
            disp('Performing linesearch, starting at s_tilde');
            lin_options = optimset('MaxIter', max_linsearch,'Display','iter','TolFun',tol_linsearch*NLAP);
            step = fminsearch(linsearch_func,0,lin_options)
        end
        s_hat = s_tilde + step.*(s_hat - s_tilde);
        beta_hat = beta_tilde + step.*(beta_hat - beta_tilde);
        NLAP_new = NLAP_func(s_hat,beta_hat);
        h_of_s_hat = h_function(s_hat);
    end
        
    objfunc_pct = (NLAP - NLAP_new)/NLAP
    s_pct = max(abs((s_tilde - s_hat)./s_tilde))
    %DEBUGGING: For debugging / crashes
    %save(['qlgeostat_iter_', num2str(num_gradevals)],'s_hat','h_of_s_hat','beta_hat','H_tilde','NLAP_new')
    
end

if NLAP_new < NLAP
    s_best = s_hat;
    beta_best = beta_hat;
    H_local = H_tilde_func(s_hat);
    NLAP = NLAP_new;    
else
    s_best = s_tilde;
    beta_best = beta_tilde;
    H_local = H_tilde;
end

finish_time = clock;
elapsed = etime(finish_time,start_time);
disp('Iteration         Time Elapsed        Obj. Func');
disp([num2str(num_gradevals), '      ', num2str(elapsed), '       ', num2str(NLAP)]);
