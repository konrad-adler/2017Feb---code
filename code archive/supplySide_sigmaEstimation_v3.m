clear all; clc;
addpath('functions/');
addpath('data/')
%cd('/Users/Lucks/Desktop/movie_project-master')
%addpath('/Users/Lucks/Desktop/movie_project-master/functions/');
%addpath('/Users/Lucks/Desktop/movie_project-master/data/');
%cd('C:\Users\Konrad\Desktop\Studium\Uni-Thesis\Movie project\2017spring - code Simon\movie_project-master - limited version')

%% Import data

% result from demand side
load ('tempdata\demandEstimation')
ntaste = Model.ntaste;
nmarket = Model.nmarket;
nmovies = Model.nmovies;

% market share data etc
load('tempdata/country_size')
load('tempdata/country_name')
load('tempdata/production_budget')
load('tempdata/release_dates','t')

% add data to model structure 
Model.country_size=country_size;
Model.budget=production_budget;

%% BACKOUT SIGMA
Model.cj_result=cj_result;
Model.mu_result=mu_result;
Model.delta_result=delta_result;
Model.gamma_result=gamma_result;
Model.sigma_result=sigma_result;

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

Model.sigma = [0.2 0.2]; % fix it for the moment

numparam=length(indexEst)*ntaste;
x0=zeros(numparam,1)+.25;
for i=1:numparam
    lb(i)= 0;
    ub(i)= 1;
end
% Estimation Loop
options = optimset('Algorithm',Model.algorithm);
options = optimset(options,'MaxIter', Model.MaxIter, 'MaxFunEvals', Model.MaxFunEvals);%,'PlotFcn',@optimplotx);
options = optimset(options,'Display', Model.MatlabDisp, 'TolFun', Model.TolFun, 'TolX', Model.TolX,'UseParallel',false);

[x_2,fval,exitflag] = fmincon(@(Params)SigmaEstimationObjective_v1(Params, Model),x0,[],[],[],[],lb,ub,...
    [],options);

[J,sigma,expprofit,likeli] = SigmaEstimationObjective_v1(x_2, Model);

bla % END OF ACTUAL CODE













