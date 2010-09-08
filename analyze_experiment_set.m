function [experiment_list experiment_data] = analyze_experiment_set( varargin )
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
ip.addRequired('input_dir',@isdir);
ip.addRequired('output_dir',@isdir);
ip.addOptional('filemask','*',@ischar);
ip.addOptional('region_marker','Position', @ischar);
ip.parse(varargin{:});

%% Loop through experiments
experiment_list = parse_experiment_dir(ip.Results.input_dir, ip.Results.filemask, ip.Results.region_marker);
experiment_data = [];

parfor e=1:size(experiment_list,2)
    experiment = experiment_list(e);
    fprintf('Analyzing Experiment: %s\n', experiment.name);
    exp_output_dir = [ip.Results.output_dir filesep experiment.name];
    [s,mess,messid] = mkdir(exp_output_dir);
    [experiment.spot_data experiment.cell_maps experiment.counts] = analyze_experiment(experiment.region_files, exp_output_dir);
    experiment_data = [experiment_data, experiment];
end

end
