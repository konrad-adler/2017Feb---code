function [J, market_share_sim] = GMMobjective(Params, Model)
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
n=Model.n;
nmarket=Model.nmarket;
nmovies=Model.nmovies;
ntaste=Model.ntaste;
W = eye(nmarket*nmovies);

%Preallocation
gammapar=zeros(ntaste,1);
c_j=zeros(nmovies*ntaste,1);
mu=zeros(nmarket*ntaste,1);
sigma=zeros(nmarket*ntaste,1);
delta_j=zeros(nmovies,1);
start_pos=1;
end_pos=ntaste;
k=1;
for i=start_pos:end_pos,
    gammapar(k)=Params(i); %Gamma guess
    k=k+1;
end

start_pos=end_pos+1;
end_pos=start_pos+nmovies*ntaste-1;
k=1;
for i=start_pos:end_pos,
    c_j(k)=Params(i); %Movie location guess
    k=k+1;
end
c_j = reshape(c_j,[nmovies,ntaste]);

start_pos=end_pos+1;
end_pos=start_pos+nmarket*ntaste-1;
k=1;
for i=start_pos:end_pos,
    mu(k)=Params(i); %Market specific consumer guess mean
    k=k+1;
end
mu = reshape(mu,[nmarket,ntaste]);
start_pos=end_pos+1;
end_pos=start_pos+nmarket*ntaste-1;
k=1;
for i=start_pos:end_pos,
    sigma(k)=Params(i); %Market specific consumer guess sigma
    k=k+1;
end
sigma = reshape(sigma,[nmarket,ntaste]);
start_pos=end_pos+1;
end_pos=start_pos+nmovies-1;
k=1;
for i=start_pos:end_pos,
    delta_j(k)=Params(i); %Delta_j
    k=k+1;
end
delta_j(1)=1;
Model.delta_j=delta_j;
%Normalizations
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
