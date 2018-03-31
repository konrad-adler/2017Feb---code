function [market_share] = Simulation_alt(Setup)
%% INPUT
n = Setup.n;
nmovies = Setup.nmovies; 
ntaste = Setup.ntaste;

market_num=Setup.j;
%X_j = Setup.X_j;
%zeta_j = Setup.zeta_j; 
c_j = Setup.c_j;


% fix coefficients
%betapar = Setup.betapar;
gammapar = Setup.gammapar;
mu = Setup.mu(market_num,:);
sigma = Setup.sigma(market_num,:);

%% PREALLOCATION

dist=zeros(n,nmovies,ntaste);
ind_utility=zeros(n,nmovies);
expenditure_share=zeros(nmovies,n);



%% CALCULATE MARKET SHARE

c_i = Setup.ctilde_i(:,market_num)*sigma + mu; % Draw consumer location

%calculate individual utility from consumer i for product j
%delta_j = X_j*betapar' + zeta_j; %Mean Utility part - identical across cons
delta_j=Setup.delta_j;

for i=1:n,
    for j=1:nmovies,
        for k=1:ntaste,
            dist(i,j,k)=gammapar(k)*abs(c_i(i,k)-c_j(j,k));
        end
    end
end

dist_sum=sum(dist,3);
for i=1:n,
    for j=1:nmovies,
        %ind_utility(i,j)= exp(dist_sum(i,j)); 
        ind_utility(i,j)= exp(delta_j(j)+dist_sum(i,j));       
    end
end

ind_utility(isinf(ind_utility))=exp(700);

ind_resist=sum(ind_utility,2);

for i=1:n,
    for j=1:nmovies,
        expenditure_share(j,i)=exp(dist_sum(i,j))/ind_resist(i);       
    end
end

market_share=mean(expenditure_share,2);

