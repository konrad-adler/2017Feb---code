function [J,Sigma_mat,expprofit,likeli] = SigmaEstimationObjective_v3(Params, Model,Results,indUtilities)
% =============================================================================================
% 
%
% INPUT: Params, vector, vector of estimated parameters
%        Model, structure
%        W, matrix, weighting matrix
% OUTPUT: J (objective)
% =============================================================================================

%% 1) DATA INPUT/PROCESS

% Fixed parameters
nmovies=Model.nmovies;
ntaste=Model.ntaste;
nmarket=Model.nmarket;
n = Model.n;

% Parameters for expected profit computations 
stepsize=Model.stepsize;
dgrid=Model.dgrid;

% Results from the demand side estimation
cj_result = Results.cj_result;
muc_result= Results.muc_result;
deltaj_result = Results.deltaj_result;
gamma_result = Results.gamma_result;
sigmac_result = Results.sigmac_result;
 
% Which movies should be estimated
indexEst = Model.indexEst;
t0 = indexEst(1)-1;

% Guess for chosen optimal taste location
cjtildeMat = reshape(Params,[length(indexEst),ntaste]); 

% Allocate variables for solution
likeli=zeros(nmovies,1);
expprofit=zeros(nmovies,1);

% variance of distances (assuming that they get it right on average and
% given the distance of all other movies)
% sigma(1)=mean(d1.^2);
% sigma(2)=mean(d2.^2);
% var-covar matrix - (assuming that there is no correlation between the
% distance in the different taste dimensions)
Sigma_mat = [Model.sigma(1) 0; 0 Model.sigma(2)];
    



%% 1) Replace the taste location of one movie taking as given the 
%   location of all other movies
guessIndUtility = ones(n,nmarket)*NaN;

for i = indexEst
    MovieIndex = i;  
    cjtilde = cjtildeMat(MovieIndex-t0,:);
    
    % what is the likelihood of observing the actual taste location given
    % the guess of the optimal taste location and the variance
    likeli(i)=mvnpdf([cj_result(MovieIndex,1) cj_result(MovieIndex,2)],cjtilde,Sigma_mat);
        
    % computing the expected profit given the guess of the optimal taste 
    % location

    % 1.1) grid around the guess
    x1 = (cjtilde(1)-dgrid):stepsize:(cjtilde(1)+dgrid);
    x2 = (cjtilde(2)-dgrid):stepsize:(cjtilde(2)+dgrid);
    gridsize=size(x1,2);
    [long,lat] = meshgrid(x1,x2);
    F = mvnpdf([long(:) lat(:)],cjtilde,Sigma_mat);
    F = reshape(F,length(long),length(lat));
    
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
            profit = Profit_OneMovie(thisIndUtility,Model,MovieIndex);%-Model.budget(i); % ENDOGENIZE BUDGET SIZE       
            tempExpProfit=tempExpProfit+ (F(i1,i2)/sum(sum(F)))*profit;
        end
    end
    expprofit(i) = tempExpProfit;
end

% save('tempdata\sigmaEst','profitMat','probaMat','tempExpProfit')
J=sum(-expprofit); 