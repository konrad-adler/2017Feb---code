function f = ind_sh_v0(expmval,expeps_ij)
% This function computes the "individual" probabilities of choosing each brand

% Written by Aviv Nevo, May 1998.

global n nmovies

% INPUT: 
% eps_ij:   n x nmovies x nmarket
% mval:     nmarket x nmovies

expmval0 = expmval';
expmval0 = repmat(expmval0,[1 1 n]);  % replicate per individual
expmval0 = permute(expmval0,[3 1 2]); % change order of dimensions

util_ij = expeps_ij.*expmval0;
util_i = sum(util_ij,2); % sum along the nmovies dimension
choiceProb = util_ij./(1+ repmat(util_i,[1 nmovies 1]));

f = choiceProb;