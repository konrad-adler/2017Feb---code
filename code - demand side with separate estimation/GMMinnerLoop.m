function [J, market_share_sim] = GMMinnerLoop(Params,mu,sigma,ctilde_i,gammapar,YearPos,t_market_share)

[nmarket, ntaste] = size(mu);
nmovies = sum(YearPos); 
market_share_sim = ones(nmovies,nmarket)*NaN;

t_c_j = Params(1:nmovies*ntaste);
t_delta_j = Params(nmovies*ntaste+1:nmovies*ntaste+nmovies);
t_c_j = reshape(t_c_j,[nmovies,ntaste]);

% normalization: first movie in the first year
if YearPos(1,1)==1
    t_delta_j(1) = 1;
end 

for j=1:nmarket,
        MarketIndex = j;
                           
        t_mu = mu(MarketIndex,:);
        t_sigma = sigma(MarketIndex,:);        
        t_ctilde_i = ctilde_i(:,MarketIndex);
        
        % simulate market shares        
        market_share_sim(:,j) = SimulationFullSample_v1(t_c_j,t_delta_j,gammapar,...
        t_mu,t_sigma,t_ctilde_i);

end

error=t_market_share-market_share_sim;
J = sum(sum(error.^2));
    



