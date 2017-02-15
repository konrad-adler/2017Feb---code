function profit = Profit_fullSample_v1(ind_utility,countrySize,MovieIndex)

[~, nmovies, ~] = size(ind_utility);

%% COMPUTE MARKET SHARES
ind_resist=sum(ind_utility,2);
ind_resist = repmat(ind_resist,[1 nmovies 1]);

expenditure_share = ind_utility./ind_resist;

market_share = mean(expenditure_share,1);
market_share = squeeze(market_share);

%% CALCULATE PROFIT
thisMovieMarketShare = market_share(MovieIndex,:);

profit = thisMovieMarketShare*countrySize';
