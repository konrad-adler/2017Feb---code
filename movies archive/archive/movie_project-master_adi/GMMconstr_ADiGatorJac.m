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

function [c,ceq] = GMMconstr_ADiGatorJac(Params,Model)
global ADiGator_GMMconstr_ADiGatorJac
if isempty(ADiGator_GMMconstr_ADiGatorJac); ADiGator_LoadData(); end
Gator1Data = ADiGator_GMMconstr_ADiGatorJac.GMMconstr_ADiGatorJac.Gator1Data;
% ADiGator Start Derivative Computations
nmarket = Model.nmarket;
%User Line: nmarket=Model.nmarket;
nmovies = Model.nmovies;
%User Line: nmovies=Model.nmovies;
ntaste = Model.ntaste;
%User Line: ntaste=Model.ntaste;
cada1f1 = nmovies*ntaste;
cada1f2 = ntaste + cada1f1;
cada1f3 = nmarket*ntaste;
pos.f = cada1f2 + cada1f3;
%User Line: pos = ntaste+nmovies*ntaste+nmarket*ntaste;
cada1f1 = pos.f + 1;
cada1f2 = nmarket*ntaste;
cada1f3 = pos.f + cada1f2;
cada1f4 = cada1f1:cada1f3;
cada1f5dx = Params.dx(Gator1Data.Index1);
cada1f5 = Params.f(cada1f4);
c.dx = -cada1f5dx;
c.f = uminus(cada1f5);
%User Line: c = -Params(pos+1:pos+nmarket*ntaste);
cada1f1 = nmovies*ntaste;
cada1f2 = ntaste + cada1f1;
cada1f3 = 2*nmarket;
cada1f4 = cada1f3*ntaste;
cada1f5 = cada1f2 + cada1f4;
cada1f6 = nmovies*ntaste;
cada1f7 = ntaste + cada1f6;
cada1f8 = 2*nmarket;
cada1f9 = cada1f8*ntaste;
cada1f10 = cada1f7 + cada1f9;
cada1f11 = cada1f10 + 3;
cada1f12 = cada1f11 + nmovies;
cada1f13 = cada1f5:cada1f12;
cada1f14dx = Params.dx(Gator1Data.Index2);
cada1f14 = Params.f(cada1f13);
cada1f15dx = -cada1f14dx;
cada1f15 = uminus(cada1f14);
cada1td1 = zeros(34,1);
cada1td1(Gator1Data.Index3) = c.dx;
cada1td1(Gator1Data.Index4) = cada1f15dx;
c.dx = cada1td1;
c.f = [c.f;cada1f15];
%User Line: c = [c;-Params(ntaste+nmovies*ntaste+2*nmarket*ntaste:ntaste+nmovies*ntaste+2*nmarket*ntaste+3+nmovies)];
ceq.f =  [];
%User Line: ceq = [];
c.dx_size = [34,64];
c.dx_location = Gator1Data.Index5;
end


function ADiGator_LoadData()
global ADiGator_GMMconstr_ADiGatorJac
ADiGator_GMMconstr_ADiGatorJac = load('GMMconstr_ADiGatorJac.mat');
return
end