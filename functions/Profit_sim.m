function profit = Profit_sim(Model,mu_result,gamma_result,cj_simul,delta,sigma_result)
% =============================================================================================
% Objective Function for Ancient city structural model
%
% INPUT: Params, vector, vector of estimated parameters
%        Model, structure
%        W, matrix, weighting matrix
% OUTPUT: J (objective)
% =============================================================================================

%% DATA INPUT/PROCESS
n=Model.n;
nmarket=Model.nmarket;
nmovies=Model.nmovies;
ntaste=Model.ntaste;
c_size=Model.country_size;
W = eye(nmarket*nmovies);
movie_replaced=Model.movie_replaced;

%Preallocation
gammapar=gamma_result;
c_j=cj_simul;
mu=mu_result;
sigma=sigma_result;
delta_j=delta;

Model.delta_j=delta_j;
%Normalizations
gammapar(1)=-gammapar(1);
gammapar(2)=-gammapar(2);

for k=1:ntaste,
    mu(Model.zerozero,k)=0; %Normalize 1 market
    mu(Model.oneone,k)=1; %Normalize 1 market
end
mu(Model.onezero,1)=0;
mu(Model.onezero,2)=1;
mu(Model.zeroone,1)=1;
mu(Model.zeroone,2)=0;


%Store in structure


% Movie unobservables/position
Model.c_j = c_j;


% fix coefficients
Model.gammapar = gammapar;
Model.mu = mu;
Model.sigma = abs(sigma);

%% AUXILIARY MODEL - SIMULATE MARKET SHARES

market_share_sim=zeros(nmovies,nmarket);
for j=1:nmarket,
    Model.j=j;
    market_share_sim(:,j)=Simulation(Model);
end

%% CALCULATE PROFIT
Profit_market=zeros(nmarket,1);
for j=1:nmarket,
    Profit_market(j)=market_share_sim(movie_replaced,j)*c_size(movie_replaced,j);
end

profit=sum(Profit_market);

end
