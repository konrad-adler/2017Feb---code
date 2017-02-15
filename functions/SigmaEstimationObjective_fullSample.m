function [J,expprofit,sigma] = SigmaEstimationObjective_fullSample(Params, Model,Results,indUtilities,c_i)
% =============================================================================================
% Given a guess for the script taste location (Params) this script computes
% - the variance between script and realized taste location
% - expected profits at the script location guess given the realized taste
% location of all movies produced previously in a given year
%
% INPUT: Params, vector, vector of estimated parameters
%        Model, structure
%        W, matrix, weighting matrix
% OUTPUT: J (objective)
% =============================================================================================

%% 1) DATA INPUT/PROCESS

% Fixed parameters
ntaste = Model.ntaste;
nmovies = Model.nmovies;
year = Model.year;
nyears = Model.nyears;
uniqueYear = unique(year);
movieNumber = 1:nmovies;

% Data
budget = Model.budget;

% Results from the demand side estimation
cj_result = Results.cj_result;
deltaj_result = Results.deltaj_result;
gammapar = Results.gamma_result;

% Guess for script location
cjtildeGuessMat = reshape(Params,[nmovies,ntaste]); 

% Allocate variables for solution
likeli=zeros(nmovies,1);
expprofit=zeros(nmovies,1);

%% 1) Compute variance between guess of script taste location and realized 
%     taste location, assuming mean zero

eps = cj_result - cjtildeGuessMat;
sigma = mean(eps.^2);

% var-covar matrix - (assuming that there is no correlation between the
% distance in the different taste dimensions)
Sigma_mat = [sigma(1) 0; 0 sigma(2)];

% note: implicit assumption is movie entrepreneus know sigma (it is fixed)


%% 2) Compute expected profits at the script location guess given the realized taste
% location of all movies produced in a given year (not sequentially for
% now)

for t = 1:nyears
    %t  =2
    YearIndex = t;
    YearPos = uniqueYear(YearIndex) == year;
    indexEst = movieNumber(YearPos);
    
    for i = 1:length(indexEst)
        %i = 4;        
        MovieIndex = indexEst(i);
        cj_obs = cj_result(MovieIndex,:);
        cjtildeGuess = cjtildeGuessMat(MovieIndex,:);
        delta_j = deltaj_result(MovieIndex); % ENDOGENIZE THIS
        
        % # ************************************************
        % IMPORTANT:
        % consider only movies in the same year
        % # ************************************************
        thisIndUtility = indUtilities(:,YearPos,:);
        
        
        tempExpProfit = UncertainProfitMap_fullSample(cjtildeGuess,Model,thisIndUtility...
            ,c_i,MovieIndex,Sigma_mat,indexEst,YearPos,gammapar,delta_j);
        
        
        expprofit(MovieIndex) = tempExpProfit - budget(MovieIndex); % ENDOGENIZE BUGDGET        
        
        % what is the likelihood of observing the actual taste location given
        % the guess of the optimal taste location and the variance
        likeli(MovieIndex)=log(mvnpdf(cj_obs,cjtildeGuess,Sigma_mat));
        
    end
end
%% 3) Return objective value

expprofit(expprofit>0)=0;
J= -sum(likeli) -sum(expprofit)/10000; 