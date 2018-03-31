function J = GMMobjective(Params, Model)
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
    zeta_j(k)=Params(i)*0; %Unobservables Movie
    k=k+1;
end


%Normalizations
gammapar(1)=-1;
gammapar(2)=-1;
mu(1,1)=5;
mu(1,2)=5;
c_j(1,1) = 1; % normalize the distance to 1
c_j(1,2) = 1; % normalize the distance to 1
zeta_j(1)=1;
sigma(1)=1;

%% AUXILIARY MODEL - SIMULATE MARKET SHARES

market_share_sim=zeros(nmovies,nmarket);
for k=1:nmarket,
    X_j = Model.X_j;
    % fix coefficients
    mu_temp = mu(k,:);
    sigma_temp = sigma(k,:);

    
    %% PREALLOCATION
    
    dist_1=zeros(n,nmovies);
    dist_2=zeros(n,nmovies);
    ind_utility=zeros(n,nmovies);
    expenditure_share=zeros(n,nmovies);
    
    
    
    %% CALCULATE MARKET SHARE
    %ctilde_i=Model.ctilde_i;
    %ctilde_i = mvnrnd(mu_temp,sigma_temp,n); % Draw consumer location
    %ctilde_i = (randn(1000,1)+mu_temp)*sigma_temp; 
    for l=1:ntaste,
        ctilde_i(:,l) = Model.ctilde_i(:,l)*sigma_temp(l) + mu_temp(l); % Draw consumer location
    end
    %ctilde_i = Model.ctilde_i(:,k)*sigma_temp + mu_temp; % Draw consumer location
    %calculate individual utility from consumer i for product j
    delta_j = X_j*betapar' + zeta_j; %Mean Utility part - identical across cons
    for i=1:n,
        for j=1:nmovies,
            dist_1(i,j)=gammapar(1)*log((ctilde_i(i,1)-c_j(j,1))^2);
        end
    end
    for i=1:n,
        for j=1:nmovies,
            dist_2(i,j)=gammapar(2)*log((ctilde_i(i,2)-c_j(j,2))^2);
        end
    end
    %dist_sum=dist_1;
    dist_sum=dist_1+dist_2;
    for i=1:n,
        for j=1:nmovies,
            ind_utility(i,j)=exp(delta_j(j)+dist_sum(i,j));
        end
    end
    
    %ind_utility(isinf(ind_utility))=exp(700);
    ind_resist=sum(ind_utility,2);
    
    for i=1:n,
        for j=1:nmovies,
            expenditure_share(i,j)=ind_utility(i,j)/ind_resist(i);
        end
    end
    market_share_sim(:,k)=sum(expenditure_share)/Model.n;
end

%% CALCULATE DISTANCE
k=1;
g=zeros(nmarket*nmovies,1);
for i=1:nmovies,
    for j=1:nmarket,
        %g(k)=abs(market_share(i,j)/market_share_sim(i,j)-1);
        g(k)=abs(market_share(i,j)-market_share_sim(i,j));
        k=k+1;
    end
end
%% OBJECTIVE FUNCTION
J = g'*W*g;


% if isnan(J)==1
%     disp(Params)
%     save('tempdata/failed','Params')
% end 
end