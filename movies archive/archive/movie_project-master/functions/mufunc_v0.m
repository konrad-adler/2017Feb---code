function f = mufunc_v0(theta2)

global n nmarket nmovies ctilde_i

% *********************************************************
% this function computes the non-linear part of the utility
%
% given a guess for: gammapar, c_j
% given draws for: ctilde_i
% *********************************************************

mu = theta2(:,1:nmarket)';
sigma = theta2(:,nmarket+1:2*nmarket)';
c_j = theta2(:,2*nmarket+1:2*nmarket+nmovies)';
gammapar = theta2(:,end);


eps_ij = ones(n,nmovies,nmarket)*NaN;

for c = 1:nmarket
    c_i = ctilde_i(:,:,c).*repmat(sigma(c,:),n,1) + repmat(mu(c,:),n,1);
    for j = 1:nmovies
        eps_ij(:,j,c) = (repmat(c_j(j,:),n,1)- c_i).^2*gammapar;      
    end
end

f = eps_ij;