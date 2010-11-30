function experiment_list = parse_experiment_list_file(filename)
% parse_experiment_list_file - creates data structure from list of
% experiments in tab delimited file:
%
% Experiment	Region	Cy3_file	Cy3.5_file	Cy5_file	DAPI_file
%
addpath bin/

experiment_list = [];
fid = fopen(filename);
data = textscan(fid, '%s\t%s\t%s\t%s\t%s\t%s', 'Delimiter', '\t');
exp_names = unique(data{1});
for e=1:size(exp_names,1)
    clear experiment
    experiment.name = exp_names{e};
    exp_rows = strcmp(data{1},experiment.name);
    experiment.regions = data{2}(exp_rows);
    experiment.region_files(:,1) = {data{3}{exp_rows}};
    experiment.region_files(:,2) = {data{4}{exp_rows}};
    experiment.region_files(:,3) = {data{5}{exp_rows}};
    experiment.region_files(:,4) = {data{6}{exp_rows}};
    experiment_list = [experiment_list experiment];
end
end