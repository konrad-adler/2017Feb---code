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
load('tempdata/balancedpanel')

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
    sigma_result(i)=x(counter);
    counter = counter+1;
end
sigma_result = reshape(sigma_result,[nmarket,ntaste]);
for i = 1:nmovies,
    fprintf( 'delta_jStar:\t %12.2f\t delta_j_guess\t %12.2f\n', 1,x(counter))
    delta_result(i)=x(counter);
    counter = counter+1;
end

bla % END OF CODE PART

save ('tempdata\demandEstimation','x','Model','gamma_result','cj_result',...
    'mu_result','sigma_result','delta_result')

%% PLOT OF THE TASTE LOCATIONS
clear all; clc;
load ('tempdata\demandEstimation')
load('tempdata/country_name')

ntaste = Model.ntaste;

figure
set(gcf,'PaperUnits','centimeters') 
xSize = 36; ySize = 36;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[100 200 xSize*50 ySize*50])
set(gcf,'Color','w')    

% normalization of countries
for k=1:ntaste,
    mu_result(Model.zerozero,k)=0; %Normalize 1 market
    mu_result(Model.oneone,k)=1; %Normalize 1 market
end
mu_result(Model.onezero,1)=0;
mu_result(Model.onezero,2)=1;
mu_result(Model.zeroone,1)=1;
mu_result(Model.zeroone,2)=0;

scatter(cj_result(:,1),cj_result(:,2),'x') % plot movies
hold on
scatter(mu_result(:,1),mu_result(:,2),'.')
hold on
for i=1:Model.nmarket,
    text(mu_result(i,1),mu_result(i,2),strrep(country_name(i),'_rev','')); % plot countries
end

print -dpng -f1 -loose output/demandEstimation.png














