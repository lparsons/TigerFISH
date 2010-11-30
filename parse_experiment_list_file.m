function experiment_set = parse_experiment_list_file(filename)
% parse_experiment_list_file - creates data structure from list of
% experiments in tab delimited file:
%
% Experiment	Region	Cy3_label   Cy3_file	Cy3.5_label Cy3.5_file	Cy5_label   Cy5_file	DAPI_label  DAPI_file
%
%       EXPERIMENT_SET - list experiment data structures
%           experiment.name - name of experiment
%           experiment.regions - list of experiment regions
%           experiment.region_files - list of files used for each region
%               (,1) = Cy3_file
%               (,2) = Cy3.5_file
%               (,3) = Cy5_file
%               (,4) = DAPI_file
%           experiment.dye_labels - cell array of dye labels for
%               each experiment
%               {Cy3_label, Cy3.5_label, Cy5_label, DAPI_label}

addpath bin/

experiment_set = [];
fid = fopen(filename);
data = textscan(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s', 'Delimiter', '\t');
exp_names = unique(data{1});
for e=1:size(exp_names,1)
    clear experiment
    experiment.name = exp_names{e};
    exp_rows = strcmp(data{1},experiment.name);
    experiment.regions = data{2}(exp_rows);
    experiment.dye_labels = [unique(data{3}(exp_rows)), unique(data{5}(exp_rows)), unique(data{7}(exp_rows)), unique(data{9}(exp_rows))];
    experiment.region_files(:,1) = {data{4}{exp_rows}};
    experiment.region_files(:,2) = {data{6}{exp_rows}};
    experiment.region_files(:,3) = {data{8}{exp_rows}};
    experiment.region_files(:,4) = {data{10}{exp_rows}};
    experiment_set = [experiment_set experiment];
end
end