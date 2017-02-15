function [J,sigma,expprofit,likeli] = SigmaEstimationObjective_v1(Params, Model)
% =============================================================================================
% 
%
% INPUT: Params, vector, vector of estimated parameters
%        Model, structure
%        W, matrix, weighting matrix
% OUTPUT: J (objective)
% =============================================================================================

%% DATA INPUT/PROCESS
nmovies=Model.nmovies;
ntaste=Model.ntaste;
nmarket = Model.nmarket;
indexEst = Model.indexEst;
t0 = indexEst(1)-1;
n = Model.n;
cj_result=Model.cj_result;

% guess for chosen optimal taste location
cjtildeMat = reshape(Params,[length(indexEst),ntaste]); 
% solution
likeli=zeros(nmovies,1);
expprofit=zeros(nmovies,1);


% variance of distances (assuming that they get it right on average and
% given the distance of all other movies)
% sigma(1)=mean(d1.^2);
% sigma(2)=mean(d2.^2);
sigma = Model.sigma;
% var-covar matrix - (assuming that there is no correlation between the
% distance in the different taste dimensions)
Sigma_mat = [sigma(1) 0; 0 sigma(2)];
    


% 1) compute individual utilities for each country/movie at realized taste
% location
indUtilities = ones(n,nmovies,nmarket)*NaN;

for j = 1:nmarket
    Model.j=j;
    for i = 1:nmovies        
        Model.movie_replaced=i;
        indUtilities(:,i,j) = Simulation_IndividualUtilities(Model);
    end
end


% 2) replace the taste location of one movie taking as given the 
%   location of all other movies
guessIndUtility = ones(n,nmarket)*NaN;

for i = indexEst
    
    Model.movie_replaced=i;
    movie_replaced=i;
    cjtilde = cjtildeMat(i-t0,:);
    
    % what is the likelihood of observing the actual taste location given
    % the guess of the optimal taste location and the variance
    likeli(i)=mvnpdf([cj_result(i,1) cj_result(i,2)],cjtilde,Sigma_mat);
        
    % computing the expected profit (for now at this location only)
    % under the guess of the optimal taste location
    Model.cj_result(movie_replaced,:) = cjtilde;
    for j = 1:nmarket
        Model.j=j;
        guessIndUtility(:,j) = Simulation_IndividualUtilities(Model);         
    end 
    
    thisIndUtility = indUtilities;
    thisIndUtility(:,movie_replaced,:) = guessIndUtility;
    
    % for profit computation consider only movies produced up to t   
    thisIndUtility = thisIndUtility(:,1:movie_replaced,:);
    
    % compute profits    
    expprofit(i) = Profit_OneMovie(thisIndUtility,Model)-Model.budget(i);    
end

expprofit=-expprofit; 
% for i=1:nmovies,
%     if expprofit(i)<0, % all movies should have positive EXPECTED payoffs
%         expprofit(i)=0;
%     end
% end
J=sum(expprofit);  %/10000+sum(-likeli); % for now not considering how likeli