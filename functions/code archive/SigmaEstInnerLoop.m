function expprofit = SigmaEstInnerLoop(Params,Model,Results,indUtilities,MovieIndex,Sigma_mat)

%% Results from the demand side estimation
muc_result= Results.muc_result;
deltaj_result = Results.deltaj_result;
gamma_result = Results.gamma_result;
sigmac_result = Results.sigmac_result;

% Fixed parameters
n = Model.n;
nmarket = Model.nmarket;

% Parameters for expected profit computations 
stepsize=Model.stepsize;
dgrid=Model.dgrid;

% Guess for optimal taste location
cjtilde = Params;

%% Grid around the guess
x1 = (cjtilde(1)-dgrid):stepsize:(cjtilde(1)+dgrid);
x2 = (cjtilde(2)-dgrid):stepsize:(cjtilde(2)+dgrid);
gridsize=size(x1,2);
[long,lat] = meshgrid(x1,x2);
F = mvnpdf([long(:) lat(:)],cjtilde,Sigma_mat);
F = reshape(F,length(long),length(lat));

%% Loop over grid
guessIndUtility = ones(n,nmarket)*NaN;
tempExpProfit=0;

for i1=1:gridsize,
    for i2=1:gridsize,
        
        c_j = [x1(i1) x2(i2)];
        delta_j = deltaj_result(MovieIndex); % ENDOGENIZE THIS
        
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
        % compute profits
        profit = Profit_OneMovie(thisIndUtility,Model,MovieIndex);
        
        % note: up to here it is the same as the profit map code where it
        % is assumed that entrepreneurs can pick a taste location freely
        tempExpProfit = tempExpProfit+ (F(i1,i2)/sum(sum(F)))*profit; 
    end    
end

expprofit = -tempExpProfit;