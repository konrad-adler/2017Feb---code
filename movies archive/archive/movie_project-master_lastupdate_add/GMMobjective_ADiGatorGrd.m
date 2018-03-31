% This code was generated using ADiGator version 1.3
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function J = GMMobjective_ADiGatorGrd(Params,Model)
global ADiGator_GMMobjective_ADiGatorGrd
if isempty(ADiGator_GMMobjective_ADiGatorGrd); ADiGator_LoadData(); end
Gator1Data = ADiGator_GMMobjective_ADiGatorGrd.GMMobjective_ADiGatorGrd.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: % =============================================================================================
%User Line: % Objective Function for Ancient city structural model
%User Line: %
%User Line: % INPUT: Params, vector, vector of estimated parameters
%User Line: %        Model, structure
%User Line: %        W, matrix, weighting matrix
%User Line: % OUTPUT: J (objective)
%User Line: % =============================================================================================
%User Line: %% DATA INPUT/PROCESS
market_share = Model.market_share;
%User Line: market_share=Model.market_share;
n = Model.n;
%User Line: n=Model.n;
nmarket = Model.nmarket;
%User Line: nmarket=Model.nmarket;
nmovies = Model.nmovies;
%User Line: nmovies=Model.nmovies;
ntaste = Model.ntaste;
%User Line: ntaste=Model.ntaste;
cada1f1 = nmarket*nmovies;
W.f = eye(cada1f1);
%User Line: W = eye(nmarket*nmovies);
%User Line: %Preallocation
gammapar.f = zeros(ntaste,1);
%User Line: gammapar=zeros(ntaste,1);
cada1f1 = nmovies*ntaste;
c_j.f = zeros(cada1f1,1);
%User Line: c_j=zeros(nmovies*ntaste,1);
cada1f1 = nmarket*ntaste;
mu.f = zeros(cada1f1,1);
%User Line: mu=zeros(nmarket*ntaste,1);
cada1f1 = nmarket*ntaste;
sigma.f = zeros(cada1f1,1);
%User Line: sigma=zeros(nmarket*ntaste,1);
betapar.f = zeros(3,1);
%User Line: betapar=zeros(3,1);
zeta_j.f = zeros(nmovies,1);
%User Line: zeta_j=zeros(nmovies,1);
start_pos.f = 1;
%User Line: start_pos=1;
posfin = ntaste;
%User Line: posfin=ntaste;
k.f = 1;
%User Line: k=1;
cadaforvar1.f = start_pos.f:posfin;
%User Line: cadaforvar1 = start_pos:posfin;
gammapar.dx = zeros(1,1);
for cadaforcount1 = 1:1
    i.f = cadaforvar1.f(:,cadaforcount1);
    %User Line: i = cadaforvar1(:,cadaforcount1);
    cada1f1dx = Params.dx(Gator1Data.Index1);
    cada1f1 = Params.f(i.f);
    gammapar.dx = cada1f1dx(Gator1Data.Index2);
    gammapar.f(k.f) = cada1f1;
    %User Line: gammapar(k)=Params(i);
    k.f = k.f + 1;
    %User Line: k=k+1;
end
start_pos.f = posfin + 1;
%User Line: start_pos=posfin+1;
cada1f1 = nmovies*ntaste;
cada1f2 = start_pos.f + cada1f1;
posfin = cada1f2 - 1;
%User Line: posfin=start_pos+nmovies*ntaste-1;
k.f = 1;
%User Line: k=1;
cadaforvar2.f = start_pos.f:posfin;
%User Line: cadaforvar2 = start_pos:posfin;
c_j.dx = zeros(20,1);
for cadaforcount2 = 1:20
    i.f = cadaforvar2.f(:,cadaforcount2);
    %User Line: i = cadaforvar2(:,cadaforcount2);
    cada1td1 = zeros(20,1);
    cada1td1(logical(Gator1Data.Index3(:,cadaforcount2))) = Params.dx(nonzeros(Gator1Data.Index3(:,cadaforcount2)));
    cada1f1dx = cada1td1;
    cada1f1 = Params.f(i.f);
    c_j.dx(logical(Gator1Data.Index4(:,cadaforcount2))) = cada1f1dx(nonzeros(Gator1Data.Index4(:,cadaforcount2)));
    c_j.f(k.f) = cada1f1;
    %User Line: c_j(k)=Params(i);
    k.f = k.f + 1;
    %User Line: k=k+1;
end
cada1f1 = [nmovies ntaste];
c_j.dx = c_j.dx;
c_j.f = reshape(c_j.f,cada1f1);
%User Line: c_j = reshape(c_j,[nmovies,ntaste]);
start_pos.f = posfin + 1;
%User Line: start_pos=posfin+1;
cada1f1 = nmarket*ntaste;
cada1f2 = start_pos.f + cada1f1;
posfin = cada1f2 - 1;
%User Line: posfin=start_pos+nmarket*ntaste-1;
k.f = 1;
%User Line: k=1;
cadaforvar3.f = start_pos.f:posfin;
%User Line: cadaforvar3 = start_pos:posfin;
mu.dx = zeros(10,1);
for cadaforcount3 = 1:10
    i.f = cadaforvar3.f(:,cadaforcount3);
    %User Line: i = cadaforvar3(:,cadaforcount3);
    cada1td1 = zeros(10,1);
    cada1td1(logical(Gator1Data.Index5(:,cadaforcount3))) = Params.dx(nonzeros(Gator1Data.Index5(:,cadaforcount3)));
    cada1f1dx = cada1td1;
    cada1f1 = Params.f(i.f);
    mu.dx(logical(Gator1Data.Index6(:,cadaforcount3))) = cada1f1dx(nonzeros(Gator1Data.Index6(:,cadaforcount3)));
    mu.f(k.f) = cada1f1;
    %User Line: mu(k)=Params(i);
    k.f = k.f + 1;
    %User Line: k=k+1;
end
cada1f1 = [nmarket ntaste];
mu.dx = mu.dx;
mu.f = reshape(mu.f,cada1f1);
%User Line: mu = reshape(mu,[nmarket,ntaste]);
start_pos.f = posfin + 1;
%User Line: start_pos=posfin+1;
cada1f1 = nmarket*ntaste;
cada1f2 = start_pos.f + cada1f1;
posfin = cada1f2 - 1;
%User Line: posfin=start_pos+nmarket*ntaste-1;
k.f = 1;
%User Line: k=1;
cadaforvar4.f = start_pos.f:posfin;
%User Line: cadaforvar4 = start_pos:posfin;
sigma.dx = zeros(10,1);
for cadaforcount4 = 1:10
    i.f = cadaforvar4.f(:,cadaforcount4);
    %User Line: i = cadaforvar4(:,cadaforcount4);
    cada1td1 = zeros(10,1);
    cada1td1(logical(Gator1Data.Index7(:,cadaforcount4))) = Params.dx(nonzeros(Gator1Data.Index7(:,cadaforcount4)));
    cada1f1dx = cada1td1;
    cada1f1 = Params.f(i.f);
    sigma.dx(logical(Gator1Data.Index8(:,cadaforcount4))) = cada1f1dx(nonzeros(Gator1Data.Index8(:,cadaforcount4)));
    sigma.f(k.f) = cada1f1;
    %User Line: sigma(k)=Params(i);
    k.f = k.f + 1;
    %User Line: k=k+1;
end
cada1f1 = [nmarket ntaste];
sigma.dx = sigma.dx;
sigma.f = reshape(sigma.f,cada1f1);
%User Line: sigma = reshape(sigma,[nmarket,ntaste]);
start_pos.f = posfin + 1;
%User Line: start_pos=posfin+1;
cada1f1 = start_pos.f + 3;
posfin = cada1f1 - 1;
%User Line: posfin=start_pos+3-1;
k.f = 1;
%User Line: k=1;
cadaforvar5.f = start_pos.f:posfin;
%User Line: cadaforvar5 = start_pos:posfin;
betapar.dx = zeros(3,1);
for cadaforcount5 = 1:3
    i.f = cadaforvar5.f(:,cadaforcount5);
    %User Line: i = cadaforvar5(:,cadaforcount5);
    cada1td1 = zeros(3,1);
    cada1td1(logical(Gator1Data.Index9(:,cadaforcount5))) = Params.dx(nonzeros(Gator1Data.Index9(:,cadaforcount5)));
    cada1f1dx = cada1td1;
    cada1f1 = Params.f(i.f);
    betapar.dx(logical(Gator1Data.Index10(:,cadaforcount5))) = cada1f1dx(nonzeros(Gator1Data.Index10(:,cadaforcount5)));
    betapar.f(k.f) = cada1f1;
    %User Line: betapar(k)=Params(i);
    k.f = k.f + 1;
    %User Line: k=k+1;
end
betapar.dx = betapar.dx;
betapar.f = betapar.f.';
%User Line: betapar=betapar';
start_pos.f = posfin + 1;
%User Line: start_pos=posfin+1;
cada1f1 = start_pos.f + nmovies;
posfin = cada1f1 - 1;
%User Line: posfin=start_pos+nmovies-1;
k.f = 1;
%User Line: k=1;
cadaforvar6.f = start_pos.f:posfin;
%User Line: cadaforvar6 = start_pos:posfin;
zeta_j.dx = zeros(20,1);
for cadaforcount6 = 1:20
    i.f = cadaforvar6.f(:,cadaforcount6);
    %User Line: i = cadaforvar6(:,cadaforcount6);
    cada1td1 = zeros(20,1);
    cada1td1(logical(Gator1Data.Index11(:,cadaforcount6))) = Params.dx(nonzeros(Gator1Data.Index11(:,cadaforcount6)));
    cada1f1dx = cada1td1;
    cada1f1 = Params.f(i.f);
    zeta_j.dx(logical(Gator1Data.Index12(:,cadaforcount6))) = cada1f1dx(nonzeros(Gator1Data.Index12(:,cadaforcount6)));
    zeta_j.f(k.f) = cada1f1;
    %User Line: zeta_j(k)=Params(i);
    k.f = k.f + 1;
    %User Line: k=k+1;
end
%User Line: %Normalizations
gammapar.f(1) = -10;
%User Line: gammapar(1)=-10;
mu.dx = mu.dx(Gator1Data.Index27);
mu.f(1) = 1;
%User Line: mu(1)=1;
c_j.dx = c_j.dx(Gator1Data.Index28);
c_j.f(1) = 1;
%User Line: c_j(1) = 1;
zeta_j.dx = zeta_j.dx(Gator1Data.Index29);
zeta_j.f(1) = 1;
%User Line: zeta_j(1)=1;
sigma.dx = sigma.dx(Gator1Data.Index30);
sigma.f(1) = 1;
%User Line: sigma(1)=1;
%User Line: %% AUXILIARY MODEL - SIMULATE MARKET SHARES
market_share_sim.f = zeros(nmovies,nmarket);
%User Line: market_share_sim=zeros(nmovies,nmarket);
cadaforvar7.f = 1:nmarket;
%User Line: cadaforvar7 = 1:nmarket;
market_share_sim.dx = zeros(8560,1);
for cadaforcount7 = 1:10
    k.f = cadaforvar7.f(:,cadaforcount7);
    %User Line: k = cadaforvar7(:,cadaforcount7);
    X_j = Model.X_j;
    %User Line: X_j = Model.X_j;
    %User Line: % fix coefficients
    cada1td1 = zeros(9,1);
    cada1td1(logical(Gator1Data.Index13(:,cadaforcount7))) = mu.dx(nonzeros(Gator1Data.Index13(:,cadaforcount7)));
    mu_temp.dx = cada1td1;
    mu_temp.f = mu.f(k.f);
    %User Line: mu_temp = mu(k);
    cada1td1 = zeros(9,1);
    cada1td1(logical(Gator1Data.Index14(:,cadaforcount7))) = sigma.dx(nonzeros(Gator1Data.Index14(:,cadaforcount7)));
    sigma_temp.dx = cada1td1;
    sigma_temp.f = sigma.f(k.f);
    %User Line: sigma_temp = sigma(k);
    %User Line: %% PREALLOCATION
    dist_1.f = zeros(50,20);
    dist_1.dx = zeros(18950,1);
    %User Line: dist_1=zeros(n,nmovies);
    ind_utility.f = zeros(50,20);
    ind_utility.dx = zeros(22500,1);
    %User Line: ind_utility=zeros(n,nmovies);
    expenditure_share.f = zeros(50,20);
    expenditure_share.dx = zeros(59000,1);
    %User Line: expenditure_share=zeros(n,nmovies);
    %User Line: %% CALCULATE MARKET SHARE
    %User Line: %ctilde_i=Model.ctilde_i;
    %User Line: %ctilde_i = mvnrnd(mu_temp,sigma_temp,n); % Draw consumer location
    %User Line: %ctilde_i = (randn(1000,1)+mu_temp)*sigma_temp;
    cada1f1 = Model.ctilde_i(:,k.f);
    cada1tempdx = sigma_temp.dx(Gator1Data.Index31);
    cada1tf1 = cada1f1(Gator1Data.Index32);
    cada1f3dx = cada1tf1(:).*cada1tempdx;
    cada1f3 = cada1f1*sigma_temp.f;
    cada1tempdx = mu_temp.dx(Gator1Data.Index33);
    cada1td1 = zeros(900,1);
    cada1td1(Gator1Data.Index34) = cada1f3dx;
    cada1td1(Gator1Data.Index35) = cada1td1(Gator1Data.Index35) + cada1tempdx;
    ctilde_i.dx = cada1td1;
    ctilde_i.f = cada1f3 + mu_temp.f;
    %User Line: ctilde_i = Model.ctilde_i(:,k)*sigma_temp + mu_temp;
    %User Line: %calculate individual utility from consumer i for product j
    cada1f1dx = betapar.dx;
    cada1f1 = betapar.f.';
    cada1f2 = 3;
    cada1td1 = zeros(3,3);
    cada1td1(Gator1Data.Index36) = cada1f1dx;
    cada1td1 = X_j*cada1td1;
    cada1td1 = cada1td1(:);
    cada1f3dx = cada1td1(Gator1Data.Index37);
    cada1f3 = X_j*cada1f1;
    cada1td1 = zeros(71,1);
    cada1td1(Gator1Data.Index38) = cada1f3dx;
    cada1td1(Gator1Data.Index39) = cada1td1(Gator1Data.Index39) + zeta_j.dx;
    delta_j.dx = cada1td1;
    delta_j.f = cada1f3 + zeta_j.f;
    %User Line: delta_j = X_j*betapar' + zeta_j;
    cadaforvar8.f = 1:n;
    %User Line: cadaforvar8 = 1:n;
    for cadaforcount8 = 1:50
        i.f = cadaforvar8.f(:,cadaforcount8);
        %User Line: i = cadaforvar8(:,cadaforcount8);
        cadaforvar9.f = 1:nmovies;
        %User Line: cadaforvar9 = 1:nmovies;
        cada1forindex1 = Gator1Data.Index15(:,cadaforcount8);
        cada1forindex5 = reshape(Gator1Data.Index22(:,cadaforcount8),18950,20);
        for cadaforcount9 = 1:20
            j.f = cadaforvar9.f(:,cadaforcount9);
            %User Line: j = cadaforvar9(:,cadaforcount9);
            cada1f1 = gammapar.f(1);
            cada1f2dx = ctilde_i.dx(cada1forindex1);
            cada1f2 = ctilde_i.f(i.f,1);
            cada1td1 = zeros(19,1);
            cada1td1(logical(Gator1Data.Index16(:,cadaforcount9))) = c_j.dx(nonzeros(Gator1Data.Index16(:,cadaforcount9)));
            cada1f3dx = cada1td1;
            cada1f3 = c_j.f(j.f,1);
            cada1td1 = zeros(37,1);
            cada1td1(Gator1Data.Index40) = cada1f2dx;
            cada1td1(Gator1Data.Index41) = cada1td1(Gator1Data.Index41) + -cada1f3dx;
            cada1f4dx = cada1td1;
            cada1f4 = cada1f2 - cada1f3;
            cada1f5dx = 2.*cada1f4.^(2-1).*cada1f4dx;
            cada1f5 = cada1f4^2;
            cada1f7dx = cada1f1.*cada1f5dx;
            cada1f7 = cada1f1*cada1f5;
            dist_1.dx(logical(cada1forindex5(:,cadaforcount9))) = cada1f7dx(nonzeros(cada1forindex5(:,cadaforcount9)));
            dist_1.f(i.f,j.f) = cada1f7;
            %User Line: dist_1(i,j)=gammapar(1)*(ctilde_i(i,1)-c_j(j,1))^2;
        end
    end
    %User Line: %     for i=1:n,
    %User Line: %         for j=1:nmovies,
    %User Line: %             dist_2(i,j)=gammapar(l)*(ctilde_i(i,2)-c_j(j,2))^2;
    %User Line: %         end
    %User Line: %     end
    dist_sum.dx = dist_1.dx;     dist_sum.f = dist_1.f;
    %User Line: dist_sum=dist_1;
    %User Line: %dist_sum=dist_1+dist_2;
    cadaforvar10.f = 1:n;
    %User Line: cadaforvar10 = 1:n;
    for cadaforcount10 = 1:50
        i.f = cadaforvar10.f(:,cadaforcount10);
        %User Line: i = cadaforvar10(:,cadaforcount10);
        cadaforvar11.f = 1:nmovies;
        %User Line: cadaforvar11 = 1:nmovies;
        cada1forindex2 = reshape(Gator1Data.Index18(:,cadaforcount10),37,20);
        cada1forindex6 = reshape(Gator1Data.Index23(:,cadaforcount10),22500,20);
        for cadaforcount11 = 1:20
            j.f = cadaforvar11.f(:,cadaforcount11);
            %User Line: j = cadaforvar11(:,cadaforcount11);
            cada1td1 = zeros(22,1);
            cada1td1(logical(Gator1Data.Index17(:,cadaforcount11))) = delta_j.dx(nonzeros(Gator1Data.Index17(:,cadaforcount11)));
            cada1f1dx = cada1td1;
            cada1f1 = delta_j.f(j.f);
            cada1td1 = zeros(37,1);
            cada1td1(logical(cada1forindex2(:,cadaforcount11))) = dist_sum.dx(nonzeros(cada1forindex2(:,cadaforcount11)));
            cada1f2dx = cada1td1;
            cada1f2 = dist_sum.f(i.f,j.f);
            cada1td1 = zeros(59,1);
            cada1td1(Gator1Data.Index42) = cada1f1dx;
            cada1td1(Gator1Data.Index43) = cada1td1(Gator1Data.Index43) + cada1f2dx;
            cada1f3dx = cada1td1;
            cada1f3 = cada1f1 + cada1f2;
            cada1f4dx = exp(cada1f3).*cada1f3dx;
            cada1f4 = exp(cada1f3);
            ind_utility.dx(logical(cada1forindex6(:,cadaforcount11))) = cada1f4dx(nonzeros(cada1forindex6(:,cadaforcount11)));
            ind_utility.f(i.f,j.f) = cada1f4;
            %User Line: ind_utility(i,j)=exp(delta_j(j)+dist_sum(i,j));
        end
    end
    %User Line: %ind_utility(isinf(ind_utility))=exp(700);
    cada1f1 = 20;
    cada1td1 = sum(sparse(Gator1Data.Index44,Gator1Data.Index45,ind_utility.dx,20,2950),1);
    ind_resist.dx = full(cada1td1(:));
    ind_resist.f = sum(ind_utility.f,2);
    %User Line: ind_resist=sum(ind_utility,2);
    cadaforvar12.f = 1:n;
    %User Line: cadaforvar12 = 1:n;
    for cadaforcount12 = 1:50
        i.f = cadaforvar12.f(:,cadaforcount12);
        %User Line: i = cadaforvar12(:,cadaforcount12);
        cadaforvar13.f = 1:nmovies;
        %User Line: cadaforvar13 = 1:nmovies;
        cada1forindex3 = reshape(Gator1Data.Index19(:,cadaforcount12),59,20);
        cada1forindex4 = Gator1Data.Index20(:,cadaforcount12);
        cada1forindex7 = reshape(Gator1Data.Index24(:,cadaforcount12),59000,20);
        for cadaforcount13 = 1:20
            j.f = cadaforvar13.f(:,cadaforcount13);
            %User Line: j = cadaforvar13(:,cadaforcount13);
            cada1td1 = zeros(59,1);
            cada1td1(logical(cada1forindex3(:,cadaforcount13))) = ind_utility.dx(nonzeros(cada1forindex3(:,cadaforcount13)));
            cada1f1dx = cada1td1;
            cada1f1 = ind_utility.f(i.f,j.f);
            cada1f2dx = ind_resist.dx(cada1forindex4);
            cada1f2 = ind_resist.f(i.f);
            cada1td1 = cada1f1dx./cada1f2;
            cada1td1 = cada1td1 + -cada1f1./cada1f2.^2.*cada1f2dx;
            cada1f3dx = cada1td1;
            cada1f3 = cada1f1/cada1f2;
            expenditure_share.dx(logical(cada1forindex7(:,cadaforcount13))) = cada1f3dx(nonzeros(cada1forindex7(:,cadaforcount13)));
            expenditure_share.f(i.f,j.f) = cada1f3;
            %User Line: expenditure_share(i,j)=ind_utility(i,j)/ind_resist(i);
        end
    end
    cada1f1 = 50;
    cada1td1 = zeros(50,1180);
    cada1td1(Gator1Data.Index46) = expenditure_share.dx;
    cada1td1 = sum(cada1td1,1);
    cada1f2dx = cada1td1(:);
    cada1f2 = sum(expenditure_share.f);
    cada1f3dx = cada1f2dx./Model.n;
    cada1f3 = cada1f2/Model.n;
    market_share_sim.dx(logical(Gator1Data.Index21(:,cadaforcount7))) = cada1f3dx(nonzeros(Gator1Data.Index21(:,cadaforcount7)));
    market_share_sim.f(:,k.f) = cada1f3;
    %User Line: market_share_sim(:,k)=sum(expenditure_share)/Model.n;
end
%User Line: %% CALCULATE DISTANCE
k.f = 1;
%User Line: k=1;
cada1f1 = nmarket*nmovies;
g.f = zeros(cada1f1,1);
%User Line: g=zeros(nmarket*nmovies,1);
cadaforvar14.f = 1:nmovies;
%User Line: cadaforvar14 = 1:nmovies;
g.dx = zeros(8560,1);
for cadaforcount14 = 1:20
    i.f = cadaforvar14.f(:,cadaforcount14);
    %User Line: i = cadaforvar14(:,cadaforcount14);
    cadaforvar15.f = 1:nmarket;
    %User Line: cadaforvar15 = 1:nmarket;
    cada1forindex1 = reshape(Gator1Data.Index25(:,cadaforcount14),59,10);
    cada1forindex2 = reshape(Gator1Data.Index26(:,cadaforcount14),8560,10);
    for cadaforcount15 = 1:10
        j.f = cadaforvar15.f(:,cadaforcount15);
        %User Line: j = cadaforvar15(:,cadaforcount15);
        cada1f1 = market_share(i.f,j.f);
        cada1td1 = zeros(59,1);
        cada1td1(logical(cada1forindex1(:,cadaforcount15))) = market_share_sim.dx(nonzeros(cada1forindex1(:,cadaforcount15)));
        cada1f2dx = cada1td1;
        cada1f2 = market_share_sim.f(i.f,j.f);
        cada1f3dx = -cada1f1./cada1f2.^2.*cada1f2dx;
        cada1f3 = cada1f1/cada1f2;
        cada1f4dx = cada1f3dx;
        cada1f4 = cada1f3 - 1;
        cada1f5dx = sign(cada1f4).*cada1f4dx;
        cada1f5 = abs(cada1f4);
        g.dx(logical(cada1forindex2(:,cadaforcount15))) = cada1f5dx(nonzeros(cada1forindex2(:,cadaforcount15)));
        g.f(k.f) = cada1f5;
        %User Line: g(k)=abs(market_share(i,j)/market_share_sim(i,j)-1);
        k.f = k.f + 1;
        %User Line: k=k+1;
    end
end
%User Line: %% OBJECTIVE FUNCTION
cada1f1dx = g.dx;
cada1f1 = g.f.';
cada1td1 = sparse(Gator1Data.Index47,Gator1Data.Index48,cada1f1dx,200,59);
cada1td1 = W.f.'*cada1td1;
cada1td1 = cada1td1(:);
cada1f2dx = full(cada1td1(Gator1Data.Index49));
cada1f2 = cada1f1*W.f;
cada1td2 = sparse(Gator1Data.Index50,Gator1Data.Index51,cada1f2dx,200,59);
cada1td2 = g.f.'*cada1td2;
cada1td1 = full(cada1td2(Gator1Data.Index52));
cada1td1 = cada1td1(:);
cada1td2 = sparse(Gator1Data.Index53,Gator1Data.Index54,g.dx,200,59);
cada1td2 = cada1f2*cada1td2;
cada1td2 = cada1td2(:);
cada1td1 = cada1td1 + full(cada1td2(Gator1Data.Index55));
J.dx = cada1td1;
J.f = cada1f2*g.f;
%User Line: J = g'*W*g;
%User Line: % if isnan(J)==1
%User Line: %     disp(Params)
%User Line: %     save('tempdata/failed','Params')
%User Line: % end
J.dx_size = 64;
J.dx_location = Gator1Data.Index56;
end


function ADiGator_LoadData()
global ADiGator_GMMobjective_ADiGatorGrd
ADiGator_GMMobjective_ADiGatorGrd = load('GMMobjective_ADiGatorGrd.mat');
return
end