function  delta_j_guess=delta_update(Model,Params)
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



%Normalizations
gammapar(1)=-1;
mu(1)=1;
c_j(1) = 0; % normalize the distance to 1

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
    market_share_sim(:,j)=Simulation_alt(Model);
end

delta_estimates=log(market_share./market_share_sim);


%% OBJECTIVE FUNCTION
delta_j_temp=mean(delta_estimates,2);
delta_j_guess=delta_j_temp/delta_j_temp(1);
delta_j_guess(1)=1;

if isnan(delta_j_guess)==1
    disp(Params)
end 
