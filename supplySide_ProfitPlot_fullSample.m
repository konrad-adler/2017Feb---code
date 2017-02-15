clear all; clc;
addpath('functions/');
addpath('data/')

%% Import data

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

%% 2) Compute profit maps
% Mapping the production into space
gridsize=50;
Model.gridsize = gridsize;
[long,lat] = meshgrid(1:1:gridsize, 1:1:gridsize);

% COMPUTE PROFIT MAP PER YEAR
year = Model.year;
nyears = Model.nyears;
uniqueYear = unique(year);

% use results from demand side estimation
cj_result = Results.cj_result;
muc_result = Results.muc_result;


for t = 1:nyears
    %t  =2
    YearIndex = t;
    YearPos = uniqueYear(YearIndex) == year;
    
    % Set the range of movies which are taken as given and the ex-ante optimal
    % location is not computed (within that year)
    temp = find(YearPos);
    t0 = temp(end)-1; % here: all movies except the last one
    yearBeg = temp(1);
    indexGiven = yearBeg:t0;
    
    MovieIndex = indexGiven(end)+1; % the next movie to be produced

    % call the function to comp
    profitMap = CertainProfitMap_fullSample(Model,Results,indUtilities,MovieIndex,indexGiven,YearPos);
    %-Model.budget(movie_replaced);
    
    smoothProfitMap=smoothn(profitMap,100);
    % find the best location for the next movie:
    bestPos = find(max(max(profitMap))==profitMap);
    [imax, jmax] = ind2sub(size(profitMap),bestPos);
    
    
    %% 3) Plot all the stuff
    figure
    set(gcf,'PaperUnits','centimeters')
    xSize = 12; ySize = 12;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[100 200 xSize*50 ySize*50])
    set(gcf,'Color','w')
    
    contour(lat,long,smoothProfitMap);
    hold on
    scatter(cj_result(indexGiven,1)*gridsize,cj_result(indexGiven,2)*gridsize) % plot movies
    for i=1:nmarket,
        scatter(muc_result(i,1)*gridsize,muc_result(i,2)*gridsize,'.'); % plot countries
        text(muc_result(i,1)*gridsize,muc_result(i,2)*gridsize,country_name(i)); % plot country names
    end
    hold on
    plot(imax,jmax,'x') % plot ideal location for next movie
    
    
    xticks = get(gca,'Xtick');
    set(gca,'xTickLabels',0:(1/(length(xticks)-1)):1)
    set(gca,'yTickLabels',0:(1/(length(xticks)-1)):1)
    title(['Profit Map for ',num2str(uniqueYear(t))])
    %bla % END OF ACTUAL CODE
    
    eval(['print -dpng -f1 -loose output/profitMap',num2str(uniqueYear(t)),'.png'])
    close all
    
end












