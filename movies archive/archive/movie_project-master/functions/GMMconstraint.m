function [c,ceq] = GMMconstraint(x,Model)
% Constraints function
ceq=x(1) + 1;
c = [];
