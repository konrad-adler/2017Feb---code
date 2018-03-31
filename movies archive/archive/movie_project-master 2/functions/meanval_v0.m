function f = meanval_v0(theta2)

% ********************************************
% This function computes the mean utility level
%
% Written by Aviv Nevo, May 1998.
% ********************************************
global s_jc

load 'tempdata/mvalold' % mvalold is exp(mvalold)

if max(abs(theta2-oldt2)) < 0.01;
	tol = 1e-9;
	flag = 0;
else
  	tol = 1e-6;
	flag = 1;
end

expeps_ij = exp(mufunc_v0(theta2));

norm = 1;
avgnorm = 1;

i = 0;
while norm > tol*10^(flag*floor(i/50)) && avgnorm > 1e-3*tol*10^(flag*floor(i/50))
    %disp(norm)
    % the following two lines are equivalent; however, the latter saves on the number of exponents
    % computed at therefore speeds up the computation by 5-10%
    %  mval = mvalold + log(s_jt) - log(mktsh(mvalold,expmu));
    mval = mvalold.*s_jc./mktsh_v0(mvalold,expeps_ij);
    %disp(mval)
    
    t = abs(mval-mvalold);
    norm = max(max(t));
    avgnorm = mean(mean(t));
    mvalold = mval;
    i = i + 1;
end
disp(['# of iterations for delta convergence:  ' num2str(i)])

if flag == 1 && max(isnan(mval)) < 1;
   mvalold = mval;
   oldt2 = theta2;
   save('tempdata/mvalold','mvalold', 'oldt2')
end   
f = log(mval);