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

% add data to model structure 
Model.country_size=country_size;
Model.budget=production_budget;

%% BACKOUT SIGMA
x_stor=x;
Model.cj_result=cj_result;
Model.mu_result=mu_result;
Model.x=x;
numparam=nmovies*2;
x0=zeros(nmovies*2,1)+.1;
clear lb ub
for i=1:numparam,
    lb(i)=-.5;
    ub(i)=.5;
end
% Estimation Loop
options = optimset('Algorithm',Model.algorithm);
options = optimset(options,'MaxIter', Model.MaxIter, 'MaxFunEvals', Model.MaxFunEvals,'PlotFcn',@optimplotx);
options = optimset(options,'Display', Model.MatlabDisp, 'TolFun', Model.TolFun, 'TolX', Model.TolX,'UseParallel',false);
[x_2,fval,exitflag] = fmincon(@(Params)GMMobjective2(Params, Model),x0,[],[],[],[],lb,ub,...
    @(Params)GMMconstr(Params, Model),options);

[J,sigma,original_loc] = GMMobjective2(x_2, Model);
%maximise likelihood, subject to expprofit positive


tic
GMMobjective2(x0, Model)
toc % 2.7 seconds

% new version
















