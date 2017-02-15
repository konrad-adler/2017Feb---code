clear all; clc;
addpath('functions/');
addpath('data/')
%cd('/Users/Lucks/Desktop/movie_project-master')
%addpath('/Users/Lucks/Desktop/movie_project-master/functions/');
%addpath('/Users/Lucks/Desktop/movie_project-master/data/');
%cd('C:\Users\Konrad\Desktop\Studium\Uni-Thesis\Movie project\2017spring - code Simon\movie_project-master - limited version')

%% 0) Import data

% result from demand side
load ('tempdata\demandEstimation')
ntaste = Model.ntaste;
nmarket = Model.nmarket;
nmovies = Model.nmovies;
n = Model.n;

% market share data etc
load('tempdata/country_size')
load('tempdata/country_name')
load('tempdata/production_budget')
load('tempdata/release_dates','t')

% add data to model structure
Model.country_size=country_size;
Model.budget=production_budget;

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

% Set the range of movies which are taken as given and the ex-ante optimal
% location is not computed
t0 = 292;
pos = find(t==t0);
indexGiven = 1:pos;
indexEst = pos+1:length(t);

Model.indexGiven = indexGiven;
Model.indexEst = indexEst;

Model.stepsize = 0.2;
Model.dgrid = 0.5;

sigmaGrid = [0.025 0.05 0.1 0.5  3];

Params = [];


for m = 1:length(sigmaGrid)
    
    Model.sigma = [sigmaGrid(m) sigmaGrid(m)]; % fix it for the moment
    
    [J,cjtildeSol,expprofit,likeli] = SigmaEstimationObjective_v4(Params, Model,Results,indUtilities);
        
    
    cjtildeMat{m} = cjtildeSol;
    likeliMat(:,m) = likeli;
    expprofitMat(:,m) = expprofit;
    
end

figure
subplot(2,1,1);plot(sum(likeliMat))
title('likelihood')
subplot(2,1,2);plot(sum(expprofitMat))
title('expected profits')

bla  % END OF CODE
%% VERSION WITHOUT LOOP BELOW










