function [market_share,mean_resist] = SimulationFullSample(Setup,MarketIndex,YearIndex)
% SIMULATE MARKET SHARES FOR ONE MARKET AND A GIVEN YEAR

%% INPUT
n = Setup.n;
ntaste = Setup.ntaste;
year = Setup.year;

YearPos = year(YearIndex) == year;
nmovies = sum(YearPos); 

c_j = Setup.c_j(YearPos,:);
delta_j=Setup.delta_j(YearPos);


% fix coefficients
gammapar = Setup.gammapar;
mu = Setup.mu(MarketIndex,:);
sigma = Setup.sigma(MarketIndex,:);

%% PREALLOCATION

dist=zeros(n,nmovies,ntaste);
ind_utility=zeros(n,nmovies);
expenditure_share=zeros(nmovies,n);
c_i=zeros(n,ntaste);


%% CALCULATE MARKET SHARE
for k=1:ntaste,
    c_i(:,k) = Setup.ctilde_i(:,MarketIndex)*sigma(k) + mu(k); % Draw consumer location
end
%calculate individual utility from consumer i for product j
%delta_j = X_j*betapar' + zeta_j; %Mean Utility part - identical across cons


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
        ind_utility(i,j)= exp(delta_j(j)+dist_sum(i,j));       
    end
end


ind_utility(isinf(ind_utility))=exp(700);
ind_resist=sum(ind_utility,2);



for i=1:n,
    for j=1:nmovies,
        expenditure_share(j,i)=ind_utility(i,j)/ind_resist(i);       
    end
end

market_share=mean(expenditure_share,2);
mean_resist=mean(ind_resist);

 

