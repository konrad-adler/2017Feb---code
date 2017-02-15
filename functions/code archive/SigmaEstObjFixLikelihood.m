function [J,cjtildeSol,expprofit] = SigmaEstObjFixLikelihood(Params, Model,Results,indUtilities)
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
muc_result= Results.muc_result;
 
% Which movies should be estimated
indexEst = Model.indexEst;
t = indexEst(1)-1;
Nest = length(indexEst);

% Allocate variables for solution
likeli=zeros(Nest,1);
expprofit=zeros(Nest,1);
cjtildeSol = zeros(Nest,ntaste);
eps = ones(Nest,ntaste)*NaN;

% variance of distances (assuming that they get it right on average and
% given the distance of all other movies)
% sigma(1)=mean(d1.^2);
% sigma(2)=mean(d2.^2);
% var-covar matrix - (assuming that there is no correlation between the
% distance in the different taste dimensions)

Sigma_mat = [Params(1) 0; 0 Params(2)];
%Sigma_mat = [Model.sigma(1) 0; 0 Model.sigma(2)];


%% 1) Replace the taste location of one movie taking as given the 
%   location of all other movies

% Options
options = optimset('Algorithm',Model.algorithm);
options = optimset(options,'MaxIter', Model.MaxIter, 'MaxFunEvals', Model.MaxFunEvals);%,'PlotFcn',@optimplotx);
options = optimset(options,'Display', Model.MatlabDisp, 'TolFun', Model.TolFun, 'TolX', Model.TolX,'UseParallel',false);


lb = [0 0];
ub = [1 1]; % TO DO: think about whether this makes sense

for i = indexEst
    %i =29;
    MovieIndex = i;  
    cj_obs = [cj_result(MovieIndex,1) cj_result(MovieIndex,2)];
    cjtilde0 = [0.5 0.5]; %cj_obs; % 
    
    % computing the expected profit given the guess of the optimal taste 
    % location
    [x,fval,~] = fmincon(@(Params)SigmaEstInnerLoop(Params,Model,Results,indUtilities,MovieIndex,Sigma_mat),...
        cjtilde0,[],[],[],[],lb,ub,@(Params)ConfidenceConstrFun(Params,cj_obs,Sigma_mat),options);
    cjtilde=x;
    tempExpProfit=fval;
    %ConfidenceConstrFun(cjtilde,cj_obs,Sigma_mat)
    
    expprofit(MovieIndex-t) = -tempExpProfit;
    cjtildeSol(MovieIndex-t,:) = cjtilde;
    eps(MovieIndex-t,:) = cj_obs - cjtilde;
    
    % what is the likelihood of observing the actual taste location given
    % the guess of the optimal taste location and the variance
    likeli(i)=log(mvnpdf(eps(MovieIndex-t,:),[0 0],Sigma_mat));
    
    % for illustration
%     figure
%     set(gcf,'Color','w')
%     scatter(cj_result(1:MovieIndex,1),cj_result(1:MovieIndex,2),'x') % plot movies
%     hold on
%     plot(cj_result(MovieIndex,1),cj_result(MovieIndex,2),'x','color','red') % plot movies
%     hold on
%     plot(cjtilde(1),cjtilde(2),'o','color','red') % plot movies    
%     hold on
%     scatter(muc_result(:,1),muc_result(:,2),'.')
%     title(['realized and optimal position ',num2str(i)])
        
end

disp('THIS ITERATION RESULTS: ')
fprintf('loglike: \t%12.2f\n',sum(likeli))
fprintf('sigma1: \t%12.2f\n',Params(1))
fprintf('sigma2: \t%12.2f\n',Params(2))
fprintf('E(eps1): \t%12.2f\n',mean(eps(:,1)))
fprintf('E(eps2): \t%12.2f\n',mean(eps(:,2)))
fprintf('var(eps1): \t%12.2f\n',var(eps(:,1)))
fprintf('var(eps2): \t%12.2f\n',var(eps(:,2)))

% save('tempdata\sigmaEst','profitMat','probaMat','tempExpProfit')

J=sum(likeli); 