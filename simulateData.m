% =========================================================================
% WRAPPER FOR SIMULATING DATA FOR TASTE SPACE ESTIMATION
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
%cd('C:\Users\Konrad\Desktop\Studium\Uni-Thesis\Movie project\2017spring - code Simon\movie_project-master')
%addpath('functions/');
cd('/Users/Lucks/Desktop/movie_project-master')
addpath('/Users/Lucks/Desktop/movie_project-master/functions/');

%% SETUP

% Global Simulation parameters
n = 100;    % individuals per market
nmarket = 6;
nmovies = 35; 
ntaste = 3; % dimension of taste space

Setup.n = n;    % individuals per market
Setup.nmarket = nmarket;
Setup.nmovies = nmovies; 
Setup.ntaste = ntaste; % dimension of taste space

% Simulate movie data


% Setup.X_j = [...
% round(unifrnd(80,120,nmovies,1)) ...% some running length
% binornd(1,0.55,nmovies,1) ...% some dummy
% exprnd(10,nmovies,1)]; % budget



% Movie unobservables/position
%Setup.zeta_j = normrnd(2,0.25,nmovies,1); 
Setup.delta_j = rand(nmovies,1); 
Setup.delta_j(1)=1;
Setup.c_j = rand(nmovies,ntaste); %Movie position in taste space
Setup.c_j(1,1) = 0; % Normalize 1 movie
Setup.c_j(1,2) = 0; % Normalize 1 movie

% fix coefficients
%Setup.betapar = [0.0005,0.0075, 0.0055];
Setup.gammapar = -rand(ntaste,1);
Setup.gammapar(1)=-1; %For now set to -1
Setup.gammapar(2)=-1; %For now set to -1
Setup.mu = rand(nmarket,ntaste);
Setup.mu(1,1)=0; %Normalize 1 market
Setup.mu(1,2)=0; %Normalize 1 market
Setup.mu(2,1)=1; %Normalize 1 market
Setup.mu(2,2)=1; %Normalize 1 market
Setup.sigma = rand(nmarket,ntaste);

% check
nmoments = nmarket*nmovies;
nparams = ntaste+nmovies*ntaste+2*nmarket*ntaste+3+nmovies;


if nmoments < nparams
    warning('Fewer moments than parameters to be estimated')
    bla
end 

% generate instruments
%???

% draw indivdual shocks for each market
ctilde_i = ones(n,nmarket)*NaN;
for i = 1:nmarket
    ctilde_i(:,i) = mvnrnd(0,1,n);
end
Setup.ctilde_i = ctilde_i;


% Simulate market specific data (position of consumers)

%% SIMULATION
market_share=zeros(nmovies,nmarket);
for j=1:nmarket,
    Setup.j=j;
    market_share(:,j)=Simulation(Setup);
end

%% OUTPUT

% save the data
save('tempdata/simulatedData','Setup','market_share')
