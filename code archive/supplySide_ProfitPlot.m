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


%% Plot country location

% Mapping the production into space
gridsize=50;
Model.gridsize = gridsize;
[long,lat] = meshgrid(1:1:gridsize, 1:1:gridsize);

% Set the range of movies which are taken as given and the ex-ante optimal
% location is not computed
t0 = 295;
pos = find(t==t0);
indexGiven = 1:pos;
Model.indexGiven = indexGiven;

% call the function to comp
profitMap = CertainProfitMap(Model,Results);

smoothProfitMap=smoothn(profitMap,100);

cj_result = Results.cj_result;
muc_result = Results.muc_result;

figure
set(gcf,'PaperUnits','centimeters') 
xSize = 12; ySize = 12;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[100 200 xSize*50 ySize*50])
set(gcf,'Color','w')    
hold on
contour(lat,long,smoothProfitMap);
hold on
scatter(cj_result(indexGiven,1)*gridsize,cj_result(indexGiven,2)*gridsize) % plot movies
for i=1:nmarket,
    text(muc_result(i,1)*gridsize,muc_result(i,2)*gridsize,country_name(i)); % plot countries
end



bla % END OF ACTUAL CODE

%%
% clear all
% load 'tempdata/workspaceNEW'
% 
% %% check difference wrt to old
% actualModel = Model;
% load ('tempdata\demandEstimation_OLD')
% clear Model
% Model = actualModel;
% 
% movie_replaced = nmovies;
% Model.movie_replaced=movie_replaced;
% 
% cjind = ntaste+1: ntaste+1+nmovies*ntaste-1;
% allcj = x(cjind);
% allcj = reshape(allcj,[nmovies,ntaste]);
% 
% allcj(movie_replaced,:) = [1/50 1/50];
% 
% x(cjind) = reshape(allcj,[length(cjind) 1]);
% Profit(x,Model)
% load 'tempdata/workspace'
% 
% 
% 
% 
% CertainProfitMap(Model,Results)
% load('tempdata\tempIndUtil')











