function [c, ceq] = GMMconstr(Params, Model)

nmarket=Model.nmarket;
nmovies=Model.nmovies;
ntaste=Model.ntaste;

pos = ntaste+nmovies*ntaste+nmarket*ntaste;

c = -Params(pos+1:pos+nmarket*ntaste);
ceq = [];

