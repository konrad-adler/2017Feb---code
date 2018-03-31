function f = mktsh_v0(expmval, expeps_ij)
% This function computes the market share for each product

% Written by Aviv Nevo, May 1998.

global n 

% ind_sh_v0 returns a matrix of dimension
% n x nmovies x nmarkets

shat_jc = sum(ind_sh_v0(expmval,expeps_ij))/n;
shat_jc = squeeze(shat_jc); % remove singleton dim

f = shat_jc';