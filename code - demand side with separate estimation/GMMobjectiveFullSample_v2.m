function [J, market_share_sim,c_j,delta_j] = GMMobjectiveFullSample_v2(Params, Model)
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
nmarket=Model.nmarket;
nmovies=Model.nmovies;
ntaste=Model.ntaste;
nyears = Model.nyears;
years = Model.year;
ctilde_i = Model.ctilde_i;
uniqueYears = unique(years);

%Preallocation
end1 = ntaste;
beg2 = end1 + 1;
end2 = beg2 + nmarket*ntaste -1;
beg3 = end2 + 1;

gammapar = Params(1:end1);
mu = Params(beg2:end2);
sigma = Params(beg3:length(Params));

mu = reshape(mu,[nmarket,ntaste]);
sigma = reshape(sigma,[nmarket,ntaste]);

% settings
options = optimset('Algorithm',Model.algorithm);
options = optimset(options,'MaxIter', Model.MaxIter, 'MaxFunEvals', Model.MaxFunEvals);
options = optimset(options,'Display', 'none', 'TolFun', Model.TolFun, 'TolX', Model.TolX,'UseParallel',false);

%% Normalizations
gammapar(1)=-gammapar(1);
gammapar(2)=-gammapar(1);

for k=1:ntaste,
    mu(Model.zerozero,k)=0; %Normalize 1 market
    mu(Model.oneone,k)=1; %Normalize 1 market
end
mu(Model.onezero,1)=0;
mu(Model.onezero,2)=1;
mu(Model.zeroone,1)=1;
mu(Model.zeroone,2)=0;


%% AUXILIARY MODEL - SIMULATE MARKET SHARES for each year

market_share_sim=zeros(nmovies,nmarket);
c_j = ones(nmovies,ntaste)*NaN;
delta_j = ones(nmovies,1)*NaN;
J = 0;

for t = 1:nyears
    YearIndex = t;
    YearPos = uniqueYears(YearIndex) == years;
    
    disp(['Start inner loop for ',num2str(uniqueYears(YearIndex))])
    % select variables for this year
    t_market_share = market_share(YearPos,:);
    t_nmovies = sum(YearPos);
    
    nparam = t_nmovies*ntaste+t_nmovies;
    InitialParamsInner = ones(nparam,1)*0.5;
    lbInner = zeros(nparam,1);
    ubInner = ones(nparam,1);
    
    [x,~,~] = fmincon(@(ParInner)GMMinnerLoop(ParInner,mu,sigma,ctilde_i,gammapar,YearPos,t_market_share)...
        ,InitialParamsInner,[],[],[],[],lbInner,ubInner,[],options);
    %x = InitialParams;
    c_j(YearPos,:) = reshape(x(1:t_nmovies*ntaste),[t_nmovies ntaste]);
    delta_j(YearPos) = x(t_nmovies*ntaste+1:end);
    
    [Jtemp, market_share_sim_temp] = GMMinnerLoop(x,mu,sigma,ctilde_i,gammapar,YearPos,t_market_share);
    
   
    market_share_sim(YearPos,:) = market_share_sim_temp;
    J = Jtemp+J;
    
end

%% OBJECTIVE FUNCTION
if isnan(J)==1
    disp(Params)
end
