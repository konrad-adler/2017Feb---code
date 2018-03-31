function f = gmmobj_v0(theta2)
% This function computes the GMM objective function

% Written by Aviv Nevo, May 1998.

global invA theta1 X_j IV nmovies nmarket

delta = meanval_v0(theta2);
delta = delta';
delta = reshape(delta,nmovies*nmarket,1); % stack markets

% the following deals with cases were the min algorithm drifts into region where the objective is not defined
if max(isnan(delta)) == 1
	f = 1e+10	   
else
	temp1 = X_j'*IV;
	temp2 = delta'*IV;
   theta1 = inv(temp1*invA*temp1')*temp1*invA*temp2';
   clear temp1 temp2 
	gmmresid = delta - X_j*theta1;
	temp1 = gmmresid'*IV;
	f = temp1*invA*temp1';
   clear temp1
   save('tempdata/gmmresid','gmmresid')
end

disp(['GMM objective:  ' num2str(f)])


