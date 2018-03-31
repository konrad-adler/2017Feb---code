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
cd('C:\Users\Konrad\Desktop\Studium\Uni-Thesis\Movie project\2017spring - code Simon\movie_project-master')
addpath('functions/');

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

% optimization settings
options = optimset('Display','iter');
options = optimset(options,'MaxIter', Model.MaxIter, 'MaxFunEvals', Model.MaxFunEvals);
options = optimset(options,'TolFun', Model.TolFun, 'TolX', Model.TolX);
%options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','iter');
%x = lsqnonlin(fun,x0,[],[],options)
%optimization
[x,fval,exitflag,output,grad] = fminunc(@(Params)GMMobjective(Params, Model),InitialParams,options);

%% compare output to 'true value'
counter = 1;

for i = 1:ntaste
    fprintf( 'gammaparStar: %12.2f gammapar %12.2f\n', gammaparStar(i),x(counter))
    counter = counter+1;
end 

for i = 1:nmovies*ntaste
    fprintf( 'c_jStar: %12.2f c_j %12.2f\n', c_jStar(i),x(counter))
    counter = counter+1;
end 

for i = 1:nmarket*ntaste
    fprintf( 'muStar: %12.2f mu %12.2f\n', muStar(i),x(counter))
    counter = counter+1;
end 

for i = 1:nmarket*ntaste
    fprintf( 'sigmaStar: %12.2f sigma %12.2f\n', sigmaStar(i),x(counter))
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



[Params, Fval, Exitflag] =  lsqnonlin(@(Params)GMMobjective(Params, Model),InitialParams,[],[],options);




