function profit = Profit_OneMovie(ind_utility,Model,MovieIndex)

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
Profit_market = zeros(nmarket,1);

for j=1:nmarket,
    Profit_market(j)=market_share(end,j)*Model.country_size(MovieIndex,j);
end


profit=sum(Profit_market);