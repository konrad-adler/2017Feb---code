function profitMap = CertainProfitMap(Model,Results,indUtilities)
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

% Fixed parameters
nmovies=Model.nmovies;
ntaste=Model.ntaste;
nmarket=Model.nmarket;
n = Model.n;
gridsize = Model.gridsize;

% Which movies should be estimated
indexGiven = Model.indexGiven;
t = indexGiven(length(indexGiven));

% Results from the demand side estimation
cj_result = Results.cj_result;
muc_result= Results.muc_result;
deltaj_result = Results.deltaj_result;
gamma_result = Results.gamma_result;
sigmac_result = Results.sigmac_result;

% Allocate variables for output
profitMap=zeros(gridsize,gridsize);


%% 1) replace the taste location of one movie taking as given the
%   location of all other movies

guessIndUtility = ones(n,nmarket)*NaN;
MovieIndex = t+1; % the next movie to be produced

for i1=1:gridsize,
    for i2=1:gridsize
        
        c_j = [i1/gridsize i2/gridsize];
        delta_j = mean(deltaj_result);%deltaj_result(MovieIndex); % ENDOGENIZE THIS
        
        for j = 1:nmarket
            MarketIndex=j;
            
            mu_c = muc_result(MarketIndex,:);
            sigma_c = sigmac_result(MarketIndex,:);
            gammapar = gamma_result;
            
            guessIndUtility(:,j) = Simulation_IndividualUtilities(Model,MarketIndex,mu_c,sigma_c,gammapar,c_j,delta_j);
        end
                
        thisIndUtility = indUtilities;
        thisIndUtility(:,MovieIndex,:) = guessIndUtility;
        
        % for profit computation consider only movies produced up to t
        thisIndUtility = thisIndUtility(:,1:MovieIndex,:);
        
        %save('tempdata\tempIndUtil','thisIndUtility')
        % compute profits
        profitMap(i1,i2) = Profit_OneMovie(thisIndUtility,Model,MovieIndex);%-Model.budget(movie_replaced);
        
    end
end


