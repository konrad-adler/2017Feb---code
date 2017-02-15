function [J,sigma,original_loc] = GMMobjective2(Params, Model)
% =============================================================================================
% Objective Function for Ancient city structural model
%
% INPUT: Params, vector, vector of estimated parameters
%        Model, structure
%        W, matrix, weighting matrix
% OUTPUT: J (objective)
% =============================================================================================

%% DATA INPUT/PROCESS
nmovies=Model.nmovies;
ntaste=Model.ntaste;
d = reshape(Params,[nmovies,ntaste]); % guess for distance to optimal taste loc
d1 = d(:,1);
d2 = d(:,2);
x=Model.x;
cj_result=Model.cj_result;
likeli=zeros(nmovies,1);
expprofit=zeros(nmovies,1);
original_loc=zeros(nmovies,2);


% variance of distances (assuming that they get it right on average and
% given the distance of all other movies)
sigma(1)=mean(d1.^2);
sigma(2)=mean(d2.^2);
% var-covar matrix - (assuming that there is no correlation between the
% distance in the different taste dimensions)
Sigma_mat = [sigma(1) 0; 0 sigma(2)];
    
    
for i=1:nmovies,
    x=Model.x;
    Model.movie_replaced=i;
    movie_replaced=i;
    % guess for optimal taste location
    mu = [cj_result(i,1)-d1(i) cj_result(i,2)-d2(i)]; 
    original_loc(i,:)=mu;
    
    
    % what is the likelihood of observing the actual taste location given
    % the guess of the optimal taste location and the variance
    likeli(i)=mvnpdf([cj_result(i,1) cj_result(i,2)],mu,Sigma_mat);
    %expprofit(i) = ExpProfit(x,Model,sigma,mu);
    
    % computing the expected profit under the guess of the optimal taste
    % location
    cjind = ntaste+1: ntaste+1+nmovies*ntaste-1;
    allcj = x(cjind);
    allcj = reshape(allcj,[nmovies,ntaste]);
    
    allcj(movie_replaced,:) = [mu(1) mu(2)];
    
    x(cjind) = reshape(allcj,[length(cjind) 1]);
    
    
%     pos=ntaste+1+(movie_replaced-1)*ntaste+1-1;
%     x(pos)=mu(1);
%     pos=ntaste+1+(movie_replaced-1)*ntaste+2-1;
%     x(pos)=mu(2);
%     
    

% this is not really expected profit - it is profit if the movie falls
    % exactly at the guess of the optimal place
    expprofit(i) = Profit(x, Model)-Model.budget(i);    
end



expprofit=-expprofit; 
for i=1:nmovies,
    if expprofit(i)<0, % all movies should have positive EXPECTED payoffs
        expprofit(i)=0;
    end
end
J=sum(expprofit)/10000+sum(-likeli);