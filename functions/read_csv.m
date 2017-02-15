function [varnames,data] = read_csv(filename,colOrder)

% this reads files where the second column is a string ("country") and all
% other nvar-2 variables are numbers

nvar = length(colOrder)/2;
fid             = fopen(filename);
varnames        = textscan(fid,'%s',nvar,'delimiter',',');



% if ind == 1
%     
% else 
%     colOrder = ['%f%s',repmat('%f',1,nvar-2)];
% end 
data            = textscan(fid,colOrder,'delimiter',',');

fclose(fid);
varnames        = varnames{1};
