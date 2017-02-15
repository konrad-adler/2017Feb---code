function [J,cjtildeSol,expprofit,sigma] = SigmaEstimationObjective_v5(Params, Model,Results,indUtilities)
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
ntaste = Model.ntaste;

% Results from the demand side estimation
cj_result = Results.cj_result;
 
% Which movies should be estimated
indexEst = Model.indexEst;
t = indexEst(1)-1;
Nest = length(indexEst);
cjtildeGuessMat = reshape(Params,[Nest,ntaste]); 

% Allocate variables for solution
likeli=zeros(Nest,1);
expprofit=zeros(Nest,1);
cjtildeSol = zeros(Nest,ntaste);

% variance of distances (assuming that they get it right on average and
% given the distance of all other movies)
disp('WRONG')
sigma(1)= mean(cjtildeGuessMat(:,1).^2);
sigma(2)= mean(cjtildeGuessMat(:,2).^2);
% var-covar matrix - (assuming that there is no correlation between the
% distance in the different taste dimensions)
Sigma_mat = [sigma(1) 0; 0 sigma(2)];


%% 1) Replace the taste location of one movie taking as given the 
%   location of all other movies

for i = indexEst
    %disp(['Movie nr: ',num2str(i)])
    %i =29;
    MovieIndex = i;  
    cj_obs = [cj_result(MovieIndex,1) cj_result(MovieIndex,2)];
    cjtildeGuess = cjtildeGuessMat(MovieIndex-t,:);
    
    tempExpProfit = SigmaEstInnerLoop(cjtildeGuess,Model,Results,indUtilities,MovieIndex,Sigma_mat);
    
    expprofit(MovieIndex-t) = (-1)*tempExpProfit - Model.budget(MovieIndex);  % inner loop returns profits*(-1),  % ENDOGENIZE BUGDGET  
    cjtildeSol(MovieIndex-t,:) = cjtildeGuess;
       
    % what is the likelihood of observing the actual taste location given
    % the guess of the optimal taste location and the variance
    likeli(MovieIndex)=log(mvnpdf(cj_obs,cjtildeGuess,Sigma_mat));
    
end

expprofit(expprofit>0)=0;


J= -sum(likeli) -sum(expprofit)/10000; 