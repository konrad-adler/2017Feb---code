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

%% Plot country location


% Profit surface: Given film and country location surface map of profits
for k=1:ntaste,
    mu_result(Model.zerozero,k)=0; %Normalize 1 market
    mu_result(Model.oneone,k)=1; %Normalize 1 market
end
mu_result(Model.onezero,1)=0;
mu_result(Model.onezero,2)=1;
mu_result(Model.zeroone,1)=1;
mu_result(Model.zeroone,2)=0;
%Which movie to replace
movie_replaced=nmovies;
Model.movie_replaced=movie_replaced;
gridsize=50;
%Mapping the production into space
[long,lat] = meshgrid(1:1:gridsize, 1:1:gridsize);

profit_map=zeros(gridsize,gridsize);
% compute profits for different movie locations - for each position in the
% taste space (movie locations are stacked as [m11 m12 m21 m22 m31 m32 ...]
% where the second digit indicates the taste space dimension)

% computation: profits at location (i1,i2) taking as given all movies
% except "movie_replaced"
% computed for the entire grid

for i1=1:gridsize,
    for i2=1:gridsize,
%         pos=ntaste+1+(movie_replaced-1)*ntaste+1-1; 
%         x(pos)=i1/50;
%         pos=ntaste+1+(movie_replaced-1)*ntaste+2-1;
%         x(pos)=i2/50;
%         
        cjind = ntaste+1: ntaste+1+nmovies*ntaste-1;
        allcj = x(cjind);
        allcj = reshape(allcj,[nmovies,ntaste]);
    
        allcj(movie_replaced,:) = [i1/50 i2/50];
    
        x(cjind) = reshape(allcj,[length(cjind) 1]);
        profit_map(i1,i2)=Profit(x,Model);
    end
    i1
end

bla1=smoothn(profit_map,100);


figure;
hold on
%surf(lat,long,bla1);
contour(lat,long,bla1);
%scatter(mu_result(:,1)*50,mu_result(:,2)*50)
hold on
scatter(cj_result(:,1)*50,cj_result(:,2)*50) % plot movies
for i=1:nmarket,
    text(mu_result(i,1)*50,mu_result(i,2)*50,country_name(i)); % plot countries
end

save('tempdata\oldprofitmap','profit_map')

bla % END OF CODE SECTION

%% PROFIT

profit_real=zeros(nmovies,1);
for i=1:nmovies,
    Model.movie_replaced=i;
    profit_real(i) = Profit(x, Model)-Model.budget(i);
end

%% EXPECTED PROFIT

%Minimal distance between movie location realisation and exp location
%How does distance map into x1 x2 sigma?

%observed location is realisation of firm entry decision such that
%exp profits are positive given a common sigma

% expected profits taking as given the position of all other movies except
% movie_replaced - NOT FOR THE ENTIRE GRID - just around mu

expprofit=zeros(nmovies,1);
for i=1:nmovies,
    mu = [cj_result(i,1) cj_result(i,2)];
    sigma(1)=.01;
    sigma(2)=.01;
    Model.movie_replaced=i;
    expprofit(i) = ExpProfit(x, Model,sigma,mu);
end


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

%% SIMULATION
beta=0.5; %add linear regression of delta on budget to get beta
%beta = regress(delta_result',log(budget))

%Draw uniformly distributed across the two dimensions sequentially, if exp
%profit larger than 0 place movie and draw location shock,
%draw until no more profitable movies

% 03/02/2017 at the moment: same as expected profit computation, but
% 1) adding a random guess of the budget which matters for mean utility
% 2) incorporating the estimated sigma

i=1;
j=1;
simnum=50;
while i<simnum,
    Model.nmovies=i;
    Model.movie_replaced=i;
    movie_replaced=i;
    mu = [rand(1) rand(1)]; % draw a random taste location for a movie
    Movie.budget(i)=rand(1)*10000000; % draw a random budget 
    delta(i)=beta*log(Movie.budget(i));  % mean utility = beta*log(budget)
    cj_simul(i,1)=mu(1);
    cj_simul(i,2)=mu(2);
    expprofit(i) = SimulProfit(Model,sigma,mu_result,gamma_result,cj_simul,delta,sigma_result,mu);
    if expprofit(i)>0
        i=i+1;
    end
    j=j+1;
    if j>80
        break
    end
end




%% Plot simulation
movie_replaced=1;
Model.movie_replaced=movie_replaced;
Model.country_size=country_size;
gridsize=50;
%Mapping the production into space
[long,lat] = meshgrid(1:1:gridsize, 1:1:gridsize);

profit_map=zeros(gridsize,gridsize);
for i1=1:gridsize,
    for i2=1:gridsize,
        cj_simul(movie_replaced,1)=i1/50;
        cj_simul(movie_replaced,2)=i2/50;
        profit_map(i1,i2)=Profit_sim(Model,mu_result,gamma_result,cj_simul,delta,sigma_result);
    end
    i1
end

bla1=smoothn(profit_map,100);
figure;
hold on
%surf(lat,long,bla1);
contour(lat,long,bla1);
hold on
scatter(cj_simul(:,1)*50,cj_simul(:,2)*50)
% Profit surface: Given film and country location surface map of profits
for k=1:ntaste,
    mu_result(Model.zerozero,k)=0; %Normalize 1 market
    mu_result(Model.oneone,k)=1; %Normalize 1 market
end
mu_result(Model.onezero,1)=0;
mu_result(Model.onezero,2)=1;
mu_result(Model.zeroone,1)=1;
mu_result(Model.zeroone,2)=0;
for i=1:nmarket,
    text(mu_result(i,1)*50,mu_result(i,2)*50,country_name(i));
end

%% COUNTERFACTUAL: EUROPEAN SUBSIDY

%Idea: Add subsidies to films that are close to Europe (roughly x1 larger
%than 25) (if mu(2)>.25 then add to expprofit


%Draw uniformly distributed across the two dimensions sequentially, if exp
%profit larger than 0 place movie and draw location shock,
%draw until no more profitable movies

i=1;
j=1;
simnum=50;
while i<simnum,
    Model.nmovies=i;
    Model.movie_replaced=i;
    movie_replaced=i;
    mu = [rand(1) rand(1)];
    Movie.budget(i)=rand(1)*10000000;
    delta(i)=beta*log(Movie.budget(i)); 
    cj_simul(i,1)=mu(1);
    cj_simul(i,2)=mu(2);
    expprofit(i) = SimulProfit(Model,sigma,mu_result,gamma_result,cj_simul,delta,sigma_result,mu);
    if mu(2)>.25
        expprofit(i)=expprofit(i)+1000000;
    end
    if expprofit(i)>0
        i=i+1;
    end
    j=j+1;
    if j>80
        break
    end
end


%% CALCULATE WELFARE NUMBERS FOR DIFFERENT COUNTRIES
n=Model.n;
nmarket=Model.nmarket;
nmovies=Model.nmovies;
ntaste=Model.ntaste;
c_size=Model.country_size;
W = eye(nmarket*nmovies);
movie_replaced=Model.movie_replaced;

%Preallocation
gammapar=gamma_result;
c_j=cj_simul;
mu=mu_result;
sigma=sigma_result;
delta_j=delta;

Model.delta_j=delta_j;
%Normalizations
gammapar(1)=-gammapar(1);
gammapar(2)=-gammapar(2);

for k=1:ntaste,
    mu(Model.zerozero,k)=0; %Normalize 1 market
    mu(Model.oneone,k)=1; %Normalize 1 market
end
mu(Model.onezero,1)=0;
mu(Model.onezero,2)=1;
mu(Model.zeroone,1)=1;
mu(Model.zeroone,2)=0;


%Store in structure


% Movie unobservables/position
Model.c_j = c_j;


% fix coefficients
Model.gammapar = gammapar;
Model.mu = mu;
Model.sigma = abs(sigma);

for j=1:nmarket,
    Model.j=j;
    welfare(j)=Simulation_Welfare(Model);
end