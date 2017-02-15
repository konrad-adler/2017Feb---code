function [gammapar, mu,delta_j] = normalization(gammapar,mu,delta_j,Model)

ntaste = Model.ntaste;

delta_j(1)=1;

gammapar(1)=-gammapar(1);
gammapar(2)=-gammapar(1);

for k=1:ntaste,
    mu(Model.zerozero,k)=0; %Normalize 1 market
    mu(Model.oneone,k)=1; %Normalize 1 market
end
mu(Model.onezero,1)=0;
mu(Model.onezero,2)=1;
mu(Model.zeroone,1)=1;
mu(Model.zeroone,2)=0;