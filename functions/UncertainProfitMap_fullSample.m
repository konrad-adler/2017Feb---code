function expprofit = UncertainProfitMap_fullSample(cjtildeGuess,Model,thisIndUtility,...
    c_i,MovieIndex,Sigma_mat,indexEst,YearPos,gammapar,delta_j)
% =============================================================================================
% COMPUTES expected profit under the assumption of all other movies
% produced in the same already exist
%
% INPUT: Params, vector, vector of estimated parameters
%        Model, structure
%        W, matrix, weighting matrix
% OUTPUT: J (objective)
% =============================================================================================

%% DATA INPUT/PROCESS

% Data
countrySize = Model.country_size(YearPos,:);
countrySize = countrySize(1,:);

% Fixed parameters
[n, ~, nmarket] = size(thisIndUtility);

stepsize=Model.stepsize;
dgrid=Model.dgrid;

t0 = indexEst(1)-1;

%% Grid around the guess
x1 = (cjtildeGuess(1)-dgrid):stepsize:(cjtildeGuess(1)+dgrid);
x2 = (cjtildeGuess(2)-dgrid):stepsize:(cjtildeGuess(2)+dgrid);
gridsize=size(x1,2);
%disp(gridsize)
[long,lat] = meshgrid(x1,x2);

%disp(size(cjtildeGuess))
F = mvnpdf([long(:) lat(:)],cjtildeGuess,Sigma_mat);
F = reshape(F,length(long),length(lat));

%% Loop over grid

guessIndUtility = ones(n,nmarket)*NaN;
tempExpProfit=0;
MovieIndexRel = MovieIndex -t0;

for i1=1:gridsize,
    for i2=1:gridsize,
        %         i1 = 1;
        %         i2 = 1;
        c_j = [x1(i1) x2(i2)];
                        
        for j = 1:nmarket
            MarketIndex=j;
            this_c_i = squeeze(c_i(:,MarketIndex,:));
            guessIndUtility(:,j) = Simulation_IndividualUtilities_Fast(gammapar,c_j,delta_j,this_c_i);
        end
        
        guessThisIndUtility = thisIndUtility;
        guessThisIndUtility(:,MovieIndexRel,:) = guessIndUtility;
                
        % compute profits        
        profit = Profit_fullSample_v1(guessThisIndUtility,countrySize,MovieIndexRel);
        
        % note: up to here it is the same as the profit map code where it
        % is assumed that entrepreneurs can pick a taste location freely
        tempExpProfit = tempExpProfit+ (F(i1,i2)/sum(sum(F)))*profit;
    end
end

expprofit = tempExpProfit;