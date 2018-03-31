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
cd('/Users/Lucks/Desktop/movie_project-master')
addpath('/Users/Lucks/Desktop/movie_project-master/functions/');

%% SETUP

% Global Simulation parameters
n = 100;    % individuals per market
nmarket = 3;
nmovies = 30; 
ntaste = 1; % dimension of taste space

Setup.n = n;    % individuals per market
Setup.nmarket = nmarket;
Setup.nmovies = nmovies; 
Setup.ntaste = ntaste; % dimension of taste space

% Simulate movie data


Setup.X_j = [...
round(unifrnd(80,120,nmovies,1)) ...% some running length
binornd(1,0.55,nmovies,1) ...% some dummy
exprnd(10,nmovies,1)]; % budget



% Movie unobservables/position
Setup.zeta_j = normrnd(2,0.25,nmovies,1); 
Setup.c_j = normrnd(0,1,nmovies,ntaste); %Movie position in taste space


% fix coefficients
Setup.betapar = [0.0005,0.0075, 0.0055];
Setup.gammapar = -rand(ntaste,1);
Setup.gammapar(1)=-1; %For now set to -1
Setup.mu = unifrnd(0,2,nmarket,ntaste);
Setup.mu(1)=1; %Normalize 1 market
Setup.sigma = unifrnd(0,2,nmarket,ntaste);
Setup.sigma(1)=1;

% generate instruments
%???


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
