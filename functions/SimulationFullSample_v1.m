function market_share = SimulationFullSample_v1(c_j,delta_j,gammapar,mu,sigma,ctilde_i)
% SIMULATE MARKET SHARES FOR ONE MARKET AND A GIVEN YEAR

%% INPUT
n = size(ctilde_i,1);
ntaste = length(gammapar);
nmovies = length(delta_j); 


%% PREALLOCATION
dist=zeros(n,nmovies,ntaste);
c_i=zeros(n,ntaste);


%% CALCULATE MARKET SHARE

for k=1:ntaste,
    c_i(:,k) = ctilde_i*sigma(k) + mu(k); % Draw consumer location
end

% calculate individual utility from consumer i for product j

for j=1:nmovies,
    for k=1:ntaste,
        dist(:,j,k)= gammapar(k)*abs(c_i(:,k)-repmat(c_j(j,k),[n,1]));
    end
end

dist_sum=sum(dist,3);

% putting together individual utility
ind_utility = exp(repmat(delta_j',[n 1]) + dist_sum);
ind_utility(isinf(ind_utility))=exp(700); % replace inf
ind_resist=sum(ind_utility,2); 

% compute individual expenditure shares
expenditure_share = ind_utility./repmat(ind_resist,[1 nmovies]);
% compute market shares
market_share=mean(expenditure_share,1)';

 

