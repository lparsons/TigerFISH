function experiment_set_data = analyze_experiment_set( varargin )
% analyze_experiment_set performs FISH image analysis to determine cell boundaries
%   and determine number of signals in each cell
%
%   analyze_experiment_set( input_dir, output_dir ) 
%       finds experiments in input_dir
%
addpath bin/

%% Parse Arguments

ip = inputParser;
ip.FunctionName = 'analyze_experiment_set';
ip.addRequired('experiment_list',@isstruct);
ip.addRequired('output_dir',@isdir);

ip.parse(varargin{:});

%% Loop through experiments
experiment_set_data = [];
parfor e=1:size(ip.Results.experiment_list,2)
    experiment = ip.Results.experiment_list(e);
    fprintf('Analyzing Experiment: %s\n', experiment.name);
    exp_output_dir = [ip.Results.output_dir filesep experiment.name];
    [s,mess,messid] = mkdir(exp_output_dir);
    [experiment.spot_data experiment.cell_maps experiment.counts] = analyze_experiment(experiment.region_files, exp_output_dir);
    experiment_set_data = [experiment_set_data, experiment];
end

%% Save results
save([ip.Results.output_dir filesep 'experiment_set_data.mat'], 'experiment_set_data')
counts = [];
for i=1:size(experiment_set_data,2)
    counts = vertcat(counts, [repmat({experiment_set_data(i).name},size(experiment_set_data(i).counts,1),1), num2cell(experiment_set_data(i).counts)]);
end
cellwrite([ip.Results.output_dir filesep 'spot_counts.csv'], counts, '\t', 'wt');
%csvwrite([ip.Results.output_dir filesep 'spot_counts.csv'], counts);
