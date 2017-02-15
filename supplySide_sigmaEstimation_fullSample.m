clear all; clc;
addpath('functions/');
addpath('data/')
%cd('/Users/Lucks/Desktop/movie_project-master')
%addpath('/Users/Lucks/Desktop/movie_project-master/functions/');
%addpath('/Users/Lucks/Desktop/movie_project-master/data/');
%cd('C:\Users\Konrad\Desktop\Studium\Uni-Thesis\Movie project\2017spring - code Simon\movie_project-master - limited version')

%% 0) Import data

% result from demand side
load ('tempdata\demandEstimation14Feb')
ntaste = Model.ntaste;
nmarket = Model.nmarket;
nmovies = Model.nmovies;
n = Model.n;

% market share data etc
load('tempdata/country_size')
load('tempdata/country_name')
load('tempdata/production_budget')
load('tempdata/release_dates')
load('tempdata/sequel')
load('tempdata/prequelT')



%% 1) compute individual utilities for each country/movie at realized taste
% location
indUtilities = ones(n,nmovies,nmarket)*NaN;

for j = 1:nmarket
    MarketIndex=j;
    for i = 1:nmovies
        MovieIndex = i;
        
        c_j = Results.cj_result(MovieIndex,:);
        mu_c = Results.muc_result(MarketIndex,:);
        sigma_c = Results.sigmac_result(MarketIndex,:);
        delta_j = Results.deltaj_result(MovieIndex);
        gammapar = Results.gamma_result;
                
        indUtilities(:,i,j) = Simulation_IndividualUtilities(Model,MarketIndex,mu_c,sigma_c,gammapar,c_j,delta_j);
    end
end


%% 2) BACKOUT SIGMA
cj_result = Results.cj_result;


%% # ************************
% 1) compute the sequel sigma
% # ************************

seqIndex = find(sequel);
preqIndex = prequelT(seqIndex); % prequelT: t of the prequel 

% drop sequels where prequel is missing
missingIndex = preqIndex==0;
seqIndex=seqIndex(~missingIndex);
preqIndex=preqIndex(~missingIndex);

% compute distance prequel - sequel taste location
nseq = length(seqIndex);
distSeqelToPrequel = ones(nseq,ntaste);

for i = 1:nseq
    distSeqelToPrequel(i,:) = cj_result(seqIndex(i),:) - cj_result(preqIndex(i),:);
end

sigmaSeq = mean(distSeqelToPrequel.^2);

%% # ***********************************************
% 2) estimate the sigma for all non-sequel movies
% # ***********************************************

% HOW DOES IT WORK:
% - we estimate the sigma from the positive expected profits condition. all movies must
% have started out from a script taste location where expected profits were
% positive
% - 1) this works as long as there are some movies for which 
% at the realized taste location profits are negative
% (force pushing for a higher sigma)
% - 2) at the same time we maximize the likelihood of observing the
% observed taste realization (force pushing for a lower sigma)

%% 2.0) drop sequels:
indUtilities = indUtilities(:,~sequel,:);
nmovies = sum(~sequel);

% add data to model structure
Model.country_size=country_size(~sequel,:);
Model.budget=production_budget(~sequel);
Model.nmovies = nmovies;
Model.year = Model.year(~sequel);

Results.cj_result = Results.cj_result(~sequel,:);
Results.deltaj_result = Results.deltaj_result(~sequel);

% compute consumer locations (to save time)
ctilde_i = Model.ctilde_i;
mu_c = Results.muc_result;
sigma_c = Results.sigmac_result;

c_i = ones(n,nmarket,ntaste)*NaN;

for k=1:ntaste
    c_i(:,:,k) = ctilde_i.*repmat(sigma_c(:,k)',[n 1]) + repmat(mu_c(:,k)',[n 1]);
end 

%% 2.1) verify the condition

% TO DO: compute profits if sigma is 0 at realized taste location and
% verify that some are < 0

%% 2.2) estimate sigma

% Parameters
Model.stepsize = 1;%0.2;
Model.dgrid = 0.1;


% Initial values
x0 = rand(nmovies*ntaste,1);%;zeros(nmovies*ntaste,1)+.25;
lb = zeros(nmovies*ntaste,1);
ub = ones(nmovies*ntaste,1);

% tic 
% [J,expprofit,sigma] = SigmaEstimationObjective_fullSample(x0, Model, Results,indUtilities,c_i);
% toc 
% bla  % END OF CODE
% Params = x0;
 

% Estimation Loop
options = optimset('Algorithm',Model.algorithm);
options = optimset(options,'MaxIter', Model.MaxIter, 'MaxFunEvals', Model.MaxFunEvals);%,'PlotFcn',@optimplotx);
options = optimset(options,'Display', Model.MatlabDisp, 'TolFun', Model.TolFun, 'TolX', Model.TolX,'UseParallel',false);

[x_2,fval,exitflag] = fmincon(@(Params)SigmaEstimationObjective_fullSample(Params, Model,Results,indUtilities,c_i),x0,[],[],[],[],lb,ub,...
    [],options);


[J,expprofit,sigma] = SigmaEstimationObjective_fullSample(x_2, Model, Results,indUtilities,c_i);


bla  % END OF CODE








