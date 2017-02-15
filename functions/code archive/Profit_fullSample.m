function profit = Profit_fullSample(ind_utility,countrySize,MovieIndex)

[n, nmovies, nmarket] = size(ind_utility);

%% COMPUTE MARKET SHARES
expenditure_share = zeros(n,nmovies,nmarket);
ind_resist=sum(ind_utility,2);


for c=1:nmarket,
    for j=1:nmovies,
        expenditure_share(:,j,c)=ind_utility(:,j,c)./ind_resist(:,1,c);
    end
end


market_share = mean(expenditure_share,1);
market_share = squeeze(market_share);

%% CALCULATE PROFIT
thisMovieMarketShare = market_share(MovieIndex,:);

profit = thisMovieMarketShare*countrySize';
