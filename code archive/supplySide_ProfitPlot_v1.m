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

% use results from demand side estimation
cj_result = Results.cj_result;
muc_result = Results.muc_result;

% find the best location for the next movie:
bestPos = find(max(max(profitMap))==profitMap);
[imax, jmax] = ind2sub(size(profitMap),bestPos);


%% plot all the stuff
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
    scatter(muc_result(i,1)*gridsize,muc_result(i,2)*gridsize,'.'); % plot countries
    text(muc_result(i,1)*gridsize,muc_result(i,2)*gridsize,country_name(i)); % plot countries
end
hold on
plot(imax,jmax,'x')

xticks = get(gca,'Xtick');
set(gca,'xTickLabels',0:(1/(length(xticks)-1)):1)
set(gca,'yTickLabels',0:(1/(length(xticks)-1)):1)

bla % END OF ACTUAL CODE

print -dpng -f1 -loose output/profitMap.png













