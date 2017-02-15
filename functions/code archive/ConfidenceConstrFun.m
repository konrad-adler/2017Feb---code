function [c, ceq] = ConfidenceConstrFun(Params,cj_obs,Sigma_mat)

% source:
% http://stats.stackexchange.com/questions/29860/confidence-interval-of-multivariate-gaussian-distribution


% Guess for optimal taste location
cjtilde = Params;

% cjtilde = cjtilde0;
% cjobs = [cj_result(MovieIndex,1) cj_result(MovieIndex,2)];

x = cjtilde';
mu = cj_obs';
k = length(x);

chistat = (x-mu)'*inv(Sigma_mat)*(x-mu);

y = chi2cdf(chistat,k);


c = y-0.95;
ceq = [];