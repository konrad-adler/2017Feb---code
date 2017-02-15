function ind_utility = Simulation_IndividualUtilities(Model,MarketIndex,mu_c,sigma_c,gammapar,c_j,delta_j)
% # **************************************
% computes individual utilities for a 
% given movie in a given coutry
% # **************************************

%% INPUT
n = Model.n;
ntaste = Model.ntaste;
ctilde_i = Model.ctilde_i;

%% PREALLOCATION
dist=zeros(n,ntaste);
ind_utility=zeros(n,1);
c_i =zeros(n,ntaste);


%% CALCULATE INDIVUDAL UTILITIES
for k=1:ntaste,
    c_i(:,k) = ctilde_i(:,MarketIndex)*sigma_c(k) + mu_c(k); % Draw consumer location
end

for i=1:n,
    for k=1:ntaste,
        dist(i,k)=gammapar(k)*abs(c_i(i,k)-c_j(1,k));
    end
end

dist_sum=sum(dist,2);

for i=1:n
    ind_utility(i)= exp(delta_j+dist_sum(i));
end

% save 'tempdata/workspaceNEW'
% bla 

% replace infinity by a very large number
ind_utility(isinf(ind_utility))=exp(700);

