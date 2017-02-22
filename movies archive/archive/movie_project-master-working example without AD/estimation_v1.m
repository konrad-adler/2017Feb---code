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
cd('/Users/Lucks/Desktop/movie_project-master')
addpath('/Users/Lucks/Desktop/movie_project-master/functions/');
%cd('C:\Users\Konrad\Desktop\Studium\Uni-Thesis\Movie project\2017spring - code Simon\movie_project-master - limited version')
%addpath('functions/');

%% IMPORT DATA
load 'tempdata/simulatedData.mat'
Model.n= Setup.n;    % individuals per market
Model.nmarket = Setup.nmarket;
Model.nmovies = Setup.nmovies; 
Model.ntaste = Setup.ntaste; % dimension of taste space
n=Model.n;
nmarket=Model.nmarket;
nmovies=Model.nmovies;
ntaste=Model.ntaste;


Model.X_j=Setup.X_j; % Import movie data
Model.market_share=market_share;
% star to save the "true" values
betaparStar = Setup.betapar;
gammaparStar = Setup.gammapar;
zeta_jStar = Setup.zeta_j;
muStar = Setup.mu;
sigmaStar = Setup.sigma;
c_jStar = Setup.c_j;

StarParams = [gammaparStar c_jStar' muStar' sigmaStar'];


%% CONFIGURATION 
Model.MaxIter=15000;                 % Optimizaiton Max Iterations
Model.MaxFunEvals=100000;             % Optimizaiton Max function evaluations
Model.TolFun=1e-14;                   % Optimizaiton Function step stopping crit
Model.TolX=1e-14;                     % Optimizaiton Control step stopping crit

%% INITIALIZATION

% initial guess
% ntaste gammaparams, nmovies*ntaste movie locations, nmarket*ntaste*2
% parameters for market specific consumer locations, 3 beta parameters
numparam=ntaste+nmovies*ntaste+2*nmarket*ntaste+3+nmovies;
%numparam = ntaste+nmovies*ntaste+2*nmarket*ntaste;

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
end_pos=start_pos+3-1;
for i=start_pos:end_pos,
    InitialParams(i)=rand(1); %Beta Params 
end
start_pos=end_pos+1;
end_pos=start_pos+nmovies-1;
for i=start_pos:end_pos,
    InitialParams(i)=rand(1); %Unobservables Movie
end

% draw indivdual shocks for each market (but only once)
ctilde_i = ones(n,nmarket)*NaN;
for i = 1:nmarket
    ctilde_i(:,i) = mvnrnd(0,1,n);
end
Model.ctilde_i = ctilde_i;

% optimization settings
options = optimset('Display','iter');
options = optimset(options,'MaxIter', Model.MaxIter, 'MaxFunEvals', Model.MaxFunEvals);
options = optimset(options,'TolFun', Model.TolFun, 'TolX', Model.TolX);
%options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','iter');
%x = lsqnonlin(fun,x0,[],[],options)

%optimization
%[x,fval,exitflag,output,grad] = fminunc(@(Params)GMMobjective(Params, Model),InitialParams,options);

options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter');%'StepTolerance',10^(-8));
[x,fval,exitflag] = fmincon(@(Params)GMMobjective(Params, Model),InitialParams,[],[],[],[],[],[],...
            @(Params)GMMconstr(Params, Model),options);

%% compare output to 'true value'
counter = 1;

for i = 1:ntaste
    if i == 1
            fprintf( 'gammaparStar:\t %12.2f\t gammapar\t normalized \n', gammaparStar(i))
    else 
        fprintf( 'gammaparStar:\t %12.2f\t gammapar\t %12.2f\n', gammaparStar(i),x(counter))
    end
    counter = counter+1;
end 

for i = 1:nmovies*ntaste
    if i == 1
        fprintf( 'c_jStar:\t %12.2f\t c_j\t  normalized \n', c_jStar(i))
    else
        fprintf( 'c_jStar:\t %12.2f\t c_j\t %12.2f\n', c_jStar(i),x(counter))
    end
    counter = counter+1;
end 

for i = 1:nmarket*ntaste
    if i == 1
        fprintf( 'muStar:\t %12.2f\t mu\t normalized \n', muStar(i))
    else
        fprintf( 'muStar:\t %12.2f\t mu\t %12.2f\n', muStar(i),x(counter))
    end
    counter = counter+1;
end 

for i = 1:nmarket*ntaste
    fprintf( 'sigmaStar:\t %12.2f\t sigma\t %12.2f\n', sigmaStar(i),x(counter))
    counter = counter+1;
end 

for i = 1:3
    fprintf( 'betaparStar: %12.2f betapar %12.2f\n', betaparStar(i),x(counter))
    counter = counter+1;
end

for i = 1:nmovies
    fprintf( 'zeta_jStar: %12.2f  zeta_j %12.2f\n', zeta_jStar(i),x(counter))
    counter = counter+1;
end 

%% some checks and various stuff
bla % end of code

[Jx, marketShareX] = GMMobjective(x, Model);
[JStar, marketShareStar] = GMMobjective(StarParams, Model);

% 1. 'sampling error': < 0.01

StarParams = [gammaparStar c_jStar' muStar' sigmaStar'];
ScaleParams = [gammaparStar c_jStar' muStar' sigmaStar'];

GMMobjective(StarParams, Model)
GMMobjective(ScaleParams, Model)



%% different optimisation 
[Params, Fval, Exitflag] =  lsqnonlin(@(Params)GMMobjective(Params, Model),InitialParams,[],[],options);




