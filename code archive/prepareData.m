clear all; clc;
addpath('functions/');
addpath('data/')


filename = 'dataEstimation_Feb2017.csv';
nvar = 58;
[varnames,data] = read_csv(filename,nvar); 

% list of all countries in the sample
country_name = {};
counter = 1;
for i = 1:length(varnames)
    if findstr('S_',varnames{i})==1
        country_name{counter} = strrep(strrep(varnames{i},'S_',''),'_rev','');
        counter = counter + 1;
    end 
    
end 
save('tempdata/country_name','country_name')    
    


for i = 1:size(data,2)
    thiscol     = data{i};
    eval([varnames{i},' = thiscol;'])
end



%% Market shares





%***************************************
% READS CSV AND SAVES THEM AS .MAT FILES
%***************************************


%% Importing Market shares
filename = 'balanced_panel.csv';
delimiter = '\t';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
balancedpanel = [dataArray{1:end-1}];
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

save ('tempdata/balancedpanel','balancedpanel')

%% Importing country size
filename = 'dataCountrySize.csv';
delimiter = '\t';
startRow = 2;
endRow = 28;
formatSpec = '%*q%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'HeaderLines', startRow-1, 'ReturnOnError', false);
fclose(fileID);
country_size = dataArray{:, 1};
clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;

save('tempdata/country_size','country_size')


filename = 'dataCountrySize.csv';
delimiter = '\t';
startRow = 2;
formatSpec = '%q%*s%[^\n\r]';


fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

country_name = dataArray{:, 1};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

save('tempdata/country_name','country_name')


%% Importing movie budget
filename = 'balanced_panel_budget.csv';
delimiter = '';
startRow = 2;
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
production_budget = dataArray{:, 1};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

save('tempdata/production_budget','production_budget')

