% =========================================================================
% WRAPPER FOR ESTIMATION ROUTINE FOR TASTE SPACE
%
% INPUT:
%
% OUTPUT:
%         Results, structure
%
% USES:
%
% =========================================================================
%% PRELIMS
clear all; clc;
%cd('/Users/Lucks/Desktop/movie_project-master')
%addpath('/Users/Lucks/Desktop/movie_project-master/functions/');
%addpath('/Users/Lucks/Desktop/movie_project-master/data/');
%cd('C:\Users\Konrad\Desktop\Studium\Uni-Thesis\Movie project\2017spring - code Simon\movie_project-master - limited version')
addpath('functions/');
addpath('data/')

%% IMPORT DATA
%load 'tempdata/simulatedData.mat'
%% Initialize variables.
filename = 'balanced_panel.csv';
delimiter = '\t';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
balancedpanel = [dataArray{1:end-1}];
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

A=corrcoef(balancedpanel);
[M,I] = min(A(:));
[I_row, I_col] = ind2sub(size(A),I);
A(I_row,I_col);
Model.zerozero=I_row; % set this as country at (0,0)
Model.oneone=I_col;   % set this as country at (1,1)
A(:,I_row)=500;
A(:,I_col)=500;
%Then find country that is both negatively correlated with min and max
% Figure out which countries are at (1,0), (0,1): which country is the
% furthest away from the two corner countries (1,1), (0,0) 
A_new=A(I_row,:)+A(I_col,:); 
[M,I] = min(A_new);
A(:,I)=500;
Model.zeroone=I;
A_new=A(I_row,:)+A(I_col,:)+A(I,:);
[M,I] = min(A_new);
Model.onezero=I;
%% Model settings
Model.n= 100;    % individuals per market
Model.nmarket = size(balancedpanel,2);
Model.nmovies = size(balancedpanel,1);
Model.ntaste = 2; % dimension of taste space
n=Model.n;
nmarket=Model.nmarket;
nmovies=Model.nmovies;
ntaste=Model.ntaste;
Model.market_share=balancedpanel;
%% CONFIGURATION
Model.MaxIter=15000;                 % Optimizaiton Max Iterations
Model.MaxFunEvals=50000;             % Optimizaiton Max function evaluations
Model.TolFun=1e-14;                   % Optimizaiton Function step stopping crit
Model.TolX=1e-14;                     % Optimizaiton Control step stopping crit
%Model.algorithm='sqp';
Model.algorithm='interior-point';
Model.MatlabDisp='iter';
%% INITIALIZATION

% draw indivdual shocks for each market (but only once)
ctilde_i = ones(n,nmarket)*NaN;
for i = 1:nmarket
    ctilde_i(:,i) = mvnrnd(0,1,n);
end
Model.ctilde_i = ctilde_i;

% InitialParams


numparam=ntaste+nmovies*ntaste+2*nmarket*ntaste+nmovies;
for i=1:numparam,
    lb(i)=0;
    ub(i)=1;
end
InitialParams=zeros(numparam,1);
start_pos=1;
end_pos=ntaste;
for i=start_pos:end_pos,
    InitialParams(i)=-rand(1); %Gamma guess
end
start_pos=end_pos+1;
end_pos=start_pos+nmovies*ntaste-1;
for i=start_pos:end_pos,
    InitialParams(i)=rand(1); %Movie location guess
end
start_pos=end_pos+1;
end_pos=start_pos+nmarket*ntaste-1;
for i=start_pos:end_pos,
    InitialParams(i)=rand(1); %Market specific consumer guess mean
end
start_pos=end_pos+1;
end_pos=start_pos+nmarket*ntaste-1;
for i=start_pos:end_pos,
    InitialParams(i)=rand(1); %Market specific consumer guess sigma
end
start_pos=end_pos+1;
end_pos=start_pos+nmovies-1;
for i=start_pos:end_pos,
    InitialParams(i)=rand(1); %Delta_j
    ub(i)=3;
end

% Estimation Loop
options = optimset('Algorithm',Model.algorithm);
options = optimset(options,'MaxIter', Model.MaxIter, 'MaxFunEvals', Model.MaxFunEvals);
options = optimset(options,'Display', Model.MatlabDisp, 'TolFun', Model.TolFun, 'TolX', Model.TolX,'UseParallel',false);
[x,fval,exitflag] = fmincon(@(Params)GMMobjective(Params, Model),InitialParams,[],[],[],[],lb,ub,...
    @(Params)GMMconstr(Params, Model),options);

bla % END OF CODE PART
save ('tempdata\demandEstimation','x','Model')

%% Output

counter = 1;

for i = 1:ntaste
    if i == 1
        fprintf( 'gammaparStar:\t %12.2f\t gammapar\t normalized \n', 1)
    else
        fprintf( 'gammaparStar:\t %12.2f\t gammapar\t %12.2f\n', 1,x(counter))
    end
    gamma_result(i)=x(counter);
    counter = counter+1;
end

for i = 1:nmovies*ntaste
    if i == 1
        fprintf( 'c_jStar:\t %12.2f\t c_j\t  normalized \n', 1)
    else
        fprintf( 'c_jStar:\t %12.2f\t c_j\t %12.2f\n', 1,x(counter))
    end
    cj_result(i)=x(counter);
    counter = counter+1;
end
cj_result = reshape(cj_result,[nmovies,ntaste]);
for i = 1:nmarket*ntaste
    if i == 1
        fprintf( 'muStar:\t %12.2f\t mu\t normalized \n', 1)
    else
        fprintf( 'muStar:\t %12.2f\t mu\t %12.2f\n', 1,x(counter))
    end
    mu_result(i)=x(counter);
    counter = counter+1;
end
mu_result = reshape(mu_result,[nmarket,ntaste]);
for i = 1:nmarket*ntaste
    fprintf( 'sigmaStar:\t %12.2f\t sigma\t %12.2f\n', 1,x(counter))
    counter = counter+1;
    sigma_result(i)=x(counter);
end
sigma_result = reshape(sigma_result,[nmarket,ntaste]);
for i = 1:nmovies,
    fprintf( 'delta_jStar:\t %12.2f\t delta_j_guess\t %12.2f\n', 1,x(counter))
    delta_result(i)=x(counter);
    counter = counter+1;
end

%% Importing country size
filename = '/Users/Lucks/Desktop/movie_project-master/data/dataCountrySize.csv';
delimiter = '\t';
startRow = 2;
endRow = 28;
formatSpec = '%*q%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'HeaderLines', startRow-1, 'ReturnOnError', false);
fclose(fileID);
country_size = dataArray{:, 1};
clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;



filename = '/Users/Lucks/Desktop/movie_project-master/data/dataCountrySize.csv';
delimiter = '\t';
startRow = 2;
formatSpec = '%q%*s%[^\n\r]';


fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);


fclose(fileID);

country_name = dataArray{:, 1};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

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
movie_replaced=1;
Model.movie_replaced=movie_replaced;
Model.country_size=country_size;
gridsize=50;
%Mapping the production into space
[long,lat] = meshgrid(1:1:gridsize, 1:1:gridsize);

profit_map=zeros(gridsize,gridsize);
for i1=1:gridsize,
    for i2=1:gridsize,
        pos=ntaste+1+(movie_replaced-1)*ntaste+1-1;
        x(pos)=i1/50;
        pos=ntaste+1+(movie_replaced-1)*ntaste+2-1;
        x(pos)=i2/50;
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
scatter(cj_result(:,1)*50,cj_result(:,2)*50)

for i=1:nmarket,
    text(mu_result(i,1)*50,mu_result(i,2)*50,country_name(i));
end


% Import movie budget
filename = '/Users/Lucks/Desktop/movie_project-master/data/balanced_panel_budget.csv';
delimiter = '';
startRow = 2;
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
production_budget = dataArray{:, 1};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;


%% PROFIT
Model.budget=production_budget;

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
