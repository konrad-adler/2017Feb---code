function [J,sigma,expprofit,likeli] = SigmaEstimationObjective_v2(Params, Model)
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


cj_result = Model.cj_result;

stepsize = Model.stepsize;
dgrid = Model.dgrid;

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
    MarketIndex=j;
    for i = 1:nmovies
        MovieIndex = i;
        
        c_j = Setup.cj_result(MovieIndex,:);
        delta_j = Setup.delta_result;
        
        % fix coefficients
        gammapar = Setup.gamma_result;
        mu = Setup.mu_result;
        sigma = Setup.sigma_result(MarketIndex,:);
        
        % normalizations
        [gammapar, mu,delta_j] = normalization(gammapar,mu,delta_j,Setup);
        
        mu = mu(MarketIndex,:);
        delta_j = delta_j(MovieIndex);
        indUtilities(:,i,j) = Simulation_IndividualUtilities_OLD(Model)
    end
end


% 2) replace the taste location of one movie taking as given the 
%   location of all other movies
guessIndUtility = ones(n,nmarket)*NaN;

for i = indexEst
    MovieIndex = i;  
    Model.movie_replaced=i;
    movie_replaced=i;
    cjtilde = cjtildeMat(i-t0,:);
    
    % what is the likelihood of observing the actual taste location given
    % the guess of the optimal taste location and the variance
    likeli(i)=mvnpdf([cj_result(i,1) cj_result(i,2)],cjtilde,Sigma_mat);
        
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
            
            Model.cj_result(movie_replaced,:) = [x1(i1) x2(i2)];            
            for j = 1:nmarket
                MarketIndex=j;
                guessIndUtility(:,j) = Simulation_IndividualUtilities(Setup,MarketIndex,mu_c,sigma_c,gammapar,c_j,delta_j)
            end            
            thisIndUtility = indUtilities;
            thisIndUtility(:,movie_replaced,:) = guessIndUtility;
            
            % for profit computation consider only movies produced up to t
            thisIndUtility = thisIndUtility(:,1:movie_replaced,:);            
            % compute profits
            profit = Profit_OneMovie(thisIndUtility,Model);%-Model.budget(i);                                    
            % sum over (density at i1,i2)/ (total density in the grid) * (profit
            % at i1,i2)
%             profitMat(i1,i2) = profit;
%             probaMat(i1,i2) = F(i1,i2)/sum(sum(F));            
            tempExpProfit=tempExpProfit+ (F(i1,i2)/sum(sum(F)))*profit;
        end
    end
    expprofit(i) = tempExpProfit;   
end

% save('tempdata\sigmaEst','profitMat','probaMat','tempExpProfit')

expprofit=-expprofit; 
% for i=1:nmovies,
%     if expprofit(i)<0, % all movies should have positive EXPECTED payoffs
%         expprofit(i)=0;
%     end
% end
J=sum(expprofit);  %/10000+sum(-likeli); % for now not considering how likeli