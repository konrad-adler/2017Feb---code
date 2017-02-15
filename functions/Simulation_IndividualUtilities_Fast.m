function ind_utility = Simulation_IndividualUtilities_Fast(gammapar,c_j,delta_j,this_c_i)
% # **************************************
% computes individual utilities for a 
% given movie in a given coutry
% # **************************************

%% CALCULATE INDIVUDAL UTILITIES
n = size(this_c_i,1);

dist_sum = abs(this_c_i-repmat(c_j,[n 1]))*gammapar';

ind_utility = exp(delta_j + dist_sum);

% replace infinity by a very large number
ind_utility(isinf(ind_utility))=exp(700);

