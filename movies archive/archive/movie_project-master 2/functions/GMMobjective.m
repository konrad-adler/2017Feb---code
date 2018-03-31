function J = GMMobjective(Params, Model)
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
betapar=zeros(3,1);
zeta_j=zeros(nmovies,1);




start_pos=1;
posfin=ntaste;
k=1;

for i=start_pos:posfin,
    gammapar(k)=Params(i); %Gamma guess
    k=k+1;
end

start_pos=posfin+1;
posfin=start_pos+nmovies*ntaste-1;
k=1;
for i=start_pos:posfin,
    c_j(k)=Params(i); %Movie location guess
    k=k+1;
end
c_j = reshape(c_j,[nmovies,ntaste]);

start_pos=posfin+1;
posfin=start_pos+nmarket*ntaste-1;
k=1;
for i=start_pos:posfin,
    mu(k)=Params(i); %Market specific consumer guess mean
    k=k+1;
end
mu = reshape(mu,[nmarket,ntaste]);
start_pos=posfin+1;
posfin=start_pos+nmarket*ntaste-1;
k=1;
for i=start_pos:posfin,
    sigma(k)=Params(i); %Market specific consumer guess sigma
    k=k+1;
end
sigma = reshape(sigma,[nmarket,ntaste]);

start_pos=posfin+1;
posfin=start_pos+3-1;
k=1;
for i=start_pos:posfin,
    betapar(k)=Params(i); %Beta Params
    k=k+1;
end
betapar=betapar';
start_pos=posfin+1;
posfin=start_pos+nmovies-1;
k=1;
for i=start_pos:posfin,
    zeta_j(k)=Params(i); %Unobservables Movie
    k=k+1;
end


%Normalizations
gammapar(1)=-1;
mu(1)=1;
sigma(1)=1;

%Store in structure


% Movie unobservables/position
Model.zeta_j = zeta_j;
Model.c_j = c_j;


% fix coefficients
Model.betapar = betapar;
Model.gammapar = gammapar;
Model.mu = mu;
Model.sigma = abs(sigma);



%% AUXILIARY MODEL - SIMULATE MARKET SHARES

market_share_sim=zeros(nmovies,nmarket);
% for k=1:nmarket,
%     Model.j=k;
%     market_share_sim(:,k)=Simulation(Model);
% end

%% CALCULATE DISTANCE

error=abs(market_share-market_share_sim);
k=1;
g=zeros(nmarket*nmovies,1);
for i=1:nmovies,
    for j=1:nmarket,
        g(k)=error(i,j);
        k=k+1;
    end
end

%% OBJECTIVE FUNCTION
J=sum(abs(g));
%object = g'*W*g;
%J=g;
end