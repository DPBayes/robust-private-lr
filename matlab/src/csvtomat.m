% Convert data from csv files to mat files
% GeneNames, GeneExpression, DrugResponse, drug_stats
% First run python hdf5tocsv.py to convert hdf5 data into csv

% Gene names
fid = fopen('../Data/ver3/GeneNames.csv');
n = textscan(fid,'%s','Delimiter', ',');
fclose(fid);
GeneNames = n{1};
save('../Data/ver3/GeneNames.mat','GeneNames');

% Gene expression
GeneExpression = csvread('../Data/ver3/GeneExpression.csv');
save('../Data/ver3/GeneExpression.mat','GeneExpression');

% Drug response
DrugResponse = csvread('../Data/ver3/DrugResponse.csv');
save('../Data/ver3/DrugResponse.mat','DrugResponse');

% Drug stats (computed with Python on a cluster)
drugstats = csvread('../Data/ver3/drugstats.csv');
sd = drugstats(1,:);
wd = drugstats(2,:);
drug_stats = struct('sd',sd,'wd',wd);
save('../wodir/drug_stats.mat','drug_stats');