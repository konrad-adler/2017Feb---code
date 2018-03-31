function [J, market_share_sim] = GMMobjectiveFullSample_v1(Params, Model)
% =============================================================================================
% Objective Function for Ancient city structural model
%
% INPUT: Params, vector, vector of estimated parameters
%        Model, structure
%        W, matrix, weighting matrix
% OUTPUT: J (objective)
% =============================================================================================

%% DATA INPUT/PROCESS

market_share=Model.market_share;
nmarket=Model.nmarket;
nmovies=Model.nmovies;
ntaste=Model.ntaste;
nyears = Model.nyears;
years = Model.year;
ctilde_i = Model.ctilde_i;
W = eye(nmarket*nmovies);
uniqueYears = unique(years);

%Preallocation
end1 = ntaste;
beg2 = end1 + 1;
end2 = beg2 + nmovies*ntaste - 1;
beg3 = end2 + 1;
end3 = beg3 + nmarket*ntaste -1; 
beg4 = end3 + 1;
end4 = beg4 + nmarket*ntaste - 1;
beg5 = end4 + 1;

gammapar = Params(1:end1);
c_j = Params(beg2:end2);
mu = Params(beg3:end3);
sigma = Params(beg4:end4);
delta_j = Params(beg5:length(Params));

c_j = reshape(c_j,[nmovies,ntaste]);
mu = reshape(mu,[nmarket,ntaste]);
sigma = reshape(sigma,[nmarket,ntaste]);

%% Normalizations
delta_j(1)=1;

gammapar(1)=-gammapar(1);
gammapar(2)=-gammapar(1);

for k=1:ntaste,
    mu(Model.zerozero,k)=0; %Normalize 1 market
    mu(Model.oneone,k)=1; %Normalize 1 market
end
mu(Model.onezero,1)=0;
mu(Model.onezero,2)=1;
mu(Model.zeroone,1)=1;
mu(Model.zeroone,2)=0;


%% AUXILIARY MODEL - SIMULATE MARKET SHARES for each year

market_share_sim=zeros(nmovies,nmarket);

for t = 1:nyears
    for j=1:nmarket,
        MarketIndex = j;
        YearIndex = t;
        YearPos = uniqueYears(YearIndex) == years;        
        
        % select parameters
        t_c_j = c_j(YearPos,:);
        t_delta_j=delta_j(YearPos);                        
        t_mu = mu(MarketIndex,:);
        t_sigma = sigma(MarketIndex,:);        
        t_ctilde_i = ctilde_i(:,MarketIndex);
        
        % simulate market shares
        market_share_sim(YearPos,j)= SimulationFullSample_v1(t_c_j,t_delta_j,gammapar,...
            t_mu,t_sigma,t_ctilde_i);
    end
end
%% CALCULATE DISTANCE
error=market_share-market_share_sim;

k=1;
g=zeros(nmarket*nmovies,1);
for i=1:nmovies,
    for j=1:nmarket,
        g(k)=error(i,j);
        k=k+1;
    end
end

%% OBJECTIVE FUNCTION
J = g'*W*g;

if isnan(J)==1
    disp(Params)
end
