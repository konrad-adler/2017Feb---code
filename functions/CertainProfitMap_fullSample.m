function profitMap = CertainProfitMap_fullSample(Model,Results,indUtilities,MovieIndex,indexGiven,YearPos)
% =============================================================================================
% COMPUTES the profit at each location in the taste space given all
% produced movies up to t0, under the assumption that a movie entrepreneur
% can exactly choose the taste location of the movie (no uncertainty)
%
% INPUT: Params, vector, vector of estimated parameters
%        Model, structure
%        W, matrix, weighting matrix
% OUTPUT: J (objective)
% =============================================================================================

%% DATA INPUT/PROCESS

% Results from the demand side estimation
muc_result= Results.muc_result;
deltaj_result = Results.deltaj_result;
gamma_result = Results.gamma_result;
sigmac_result = Results.sigmac_result;

% Data
countrySize = Model.country_size(YearPos,:);
countrySize = countrySize(1,:);

% Fixed parameters
nmarket=Model.nmarket;
n = Model.n;
gridsize = Model.gridsize;

t0 = indexGiven(1)-1;

% Allocate variables for output
profitMap=zeros(gridsize,gridsize);


%% Loop over grid
delta_j = mean(deltaj_result(YearPos));%deltaj_result(MovieIndex); % ENDOGENIZE THIS

guessIndUtility = ones(n,nmarket)*NaN;


for i1=1:gridsize,
    for i2=1:gridsize
        
        c_j = [i1/gridsize i2/gridsize];        
        
        for j = 1:nmarket
            MarketIndex=j;
            
            mu_c = muc_result(MarketIndex,:);
            sigma_c = sigmac_result(MarketIndex,:);
            gammapar = gamma_result;
            
            guessIndUtility(:,j) = Simulation_IndividualUtilities(Model,MarketIndex,mu_c,sigma_c,gammapar,c_j,delta_j);
        end
                
        thisIndUtility = indUtilities;
        thisIndUtility(:,MovieIndex,:) = guessIndUtility;
        
        % # ************************************************
        % IMPORTANT:
        % for profit computation consider movies that have already been
        % produced this year (indexGiven) + the next one
        % # ************************************************
        thisIndUtility = thisIndUtility(:,[indexGiven MovieIndex],:);                
        % compute profits
        MovieIndexRel = MovieIndex -t0;
        profitMap(i1,i2) = Profit_fullSample(thisIndUtility,countrySize,MovieIndexRel);
        
    end
end


