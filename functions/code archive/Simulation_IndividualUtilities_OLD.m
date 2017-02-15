function ind_utility = Simulation_IndividualUtilities_OLD(Setup,MovieIndex,MarketIndex)
% # **************************************
%% computes individual utilities for a given movie in a given coutry

%% INPUT
n = Setup.n;
ntaste = Setup.ntaste;

% this movie and country
market_num = MarketIndex;
movie_num  = MovieIndex;

c_j = Setup.cj_result(movie_num,:);
delta_j = Setup.delta_result;

% fix coefficients
gammapar = Setup.gamma_result;
mu = Setup.mu_result;
sigma = Setup.sigma_result(market_num,:);

% normalizations
[gammapar, mu,delta_j] = normalization(gammapar,mu,delta_j,Setup);

mu = mu(market_num,:);
delta_j = delta_j(movie_num);
%% PREALLOCATION
dist=zeros(n,ntaste);
ind_utility=zeros(n,1);
c_i =zeros(n,ntaste);


%% CALCULATE INDIVUDAL UTILITIES
for k=1:ntaste,
    c_i(:,k) = Setup.ctilde_i(:,market_num)*sigma(k) + mu(k); % Draw consumer location
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

% replace infinity by a very large number
ind_utility(isinf(ind_utility))=exp(700);

%save 'tempdata/workspace_new'
