clear all; clc;
addpath('functions/');
addpath('data/')

%% read the csv

filename = 'dataEstimation_Feb2017_full.csv';
colOrder = [repmat('%f',1,58) repmat('%s',1,2) repmat('%f',1,1)];

%colOrder = [repmat('%f',1,58)];


[varnames,data] = read_csv(filename,colOrder); 


%% list of all countries in the sample
country_name = {};
counter = 1;
for i = 1:length(varnames)
    if findstr('S_',varnames{i})==1
        country_name{counter} = strrep(strrep(varnames{i},'S_',''),'_rev','');
        counter = counter + 1;
    end 
    
end 
save('tempdata/country_name','country_name')    

nmarket = length(country_name);

%% create variables

for i = 1:size(data,2)
    thiscol     = data{i};
    eval([varnames{i},' = thiscol;'])
end

%% market shares 

balancedpanel = ones(length(t),nmarket)*NaN;

for i = 1:nmarket
    eval(['balancedpanel(:,i) = S_',country_name{i},'_rev;'])
end 

save ('tempdata/balancedpanel','balancedpanel')

%% country size

country_size = ones(length(t),nmarket)*NaN;

for i = 1:nmarket
    eval(['country_size(:,i) = X_',country_name{i},'_rev;'])
end 

save ('tempdata/country_size','country_size')

%% bunch of other variables

% release dates
save('tempdata/release_dates','t')
% budget size
save('tempdata/production_budget','production_budget')
% sequel
save('tempdata/sequel','sequel')
% year
save('tempdata/year','year')
% 
save('tempdata/title','title')
save('tempdata/prequelT','prequelT')

