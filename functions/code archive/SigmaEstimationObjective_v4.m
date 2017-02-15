function [J,cjtildeSol,expprofit,likeli] = SigmaEstimationObjective_v4(Params, Model,Results,indUtilities)
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
ntaste = Model.ntaste;

% Results from the demand side estimation
cj_result = Results.cj_result;
 
% Which movies should be estimated
indexEst = Model.indexEst;
t = indexEst(1)-1;
% Allocate variables for solution
likeli=zeros(nmovies,1);
expprofit=zeros(nmovies,1);
cjtildeSol = zeros(length(indexEst),ntaste);

% variance of distances (assuming that they get it right on average and
% given the distance of all other movies)
% sigma(1)=mean(d1.^2);
% sigma(2)=mean(d2.^2);
% var-covar matrix - (assuming that there is no correlation between the
% distance in the different taste dimensions)
Sigma_mat = [Model.sigma(1) 0; 0 Model.sigma(2)];


%% 1) Replace the taste location of one movie taking as given the 
%   location of all other movies

% Options
options = optimset('Algorithm',Model.algorithm);
options = optimset(options,'MaxIter', Model.MaxIter, 'MaxFunEvals', Model.MaxFunEvals);%,'PlotFcn',@optimplotx);
options = optimset(options,'Display', Model.MatlabDisp, 'TolFun', Model.TolFun, 'TolX', Model.TolX,'UseParallel',false);

cjtilde0 = [0.25 0.25];
lb = [0 0];
ub = [1 1];

for i = indexEst
    MovieIndex = i;  
                
    % computing the expected profit given the guess of the optimal taste 
    % location
    [cjtilde,tempExpProfit,~] = fmincon(@(Params)SigmaEstInnerLoop(Params,Model,Results,indUtilities,MovieIndex,Sigma_mat),...
        cjtilde0,[],[],[],[],lb,ub,[],options);
    
    expprofit(MovieIndex) = -tempExpProfit;
    cjtildeSol(MovieIndex-t,:) = cjtilde;
    
    % what is the likelihood of observing the actual taste location given
    % the guess of the optimal taste location and the variance
    likeli(i)=mvnpdf([cj_result(MovieIndex,1) cj_result(MovieIndex,2)],cjtilde,Sigma_mat);
        
end

% save('tempdata\sigmaEst','profitMat','probaMat','tempExpProfit')
J=sum(expprofit); 