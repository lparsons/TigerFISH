function experiment_counts = analyze_experiment_set( experiment_list, output_dir, varargin )
% analyze_experiment_set performs FISH image analysis to determine cell boundaries
%   and determine number of signals in each cell
%
%   experiment_set_data = analyze_experiment_set( experiment_list, output_dir, 'ParamName', ParamValue, ... ) 
%
%   INPUT
%       experiment_list - list of experiment data structures
%           experiment.name - name of experiment
%           experiment.regions - list of experiment regions
%           experiment.region_files - list of files used for each region
%               (,1) = Cy3_file
%               (,2) = Cy3.5_file
%               (,3) = Cy5_file
%               (,4) = DAPI_file	
%
%       output_dir - Output directory used to save results, intermediate
%           data, etc.
%
%
%   OPTIONAL PARAMETERS
%       params - optional struct containing parameter values for analysis
%           See default_parameters.m for list of parameters and
%           documentation
%
%       load_results - optional parameter, if true load previous cell map and
%           spot intensity data (if it exists).  Off be default
%
%
%   OUTPUT
%       experiment_set_data - struct with the following fields
%
%           name - experiment name
%           regions - list of regions
%           region_files - for each region, list of cy3, cy3.5, cy5, dapi
%               files
%           spot_data - struct with fields for each cy dye
%               each dye contains a 7 column matrix with row for each spot
%               columns are: region, x, y, z, intensity, cell, cell_type
%                   cell_type 0=normal, 1=edge_cell, 2=background (not cell)
%           cell_maps - cell array with 2D matrix identifying cell pixels
%           counts - matrix with 5 colums, row for each cell
%               columns: region, cell, cy3 count, cy3.5 count, cy5 count
%
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------
p = mfilename('fullpath');
[pathstr] = fileparts(p);
addpath([pathstr filesep 'bin'])

%% Parse Arguments
ip = inputParser;
ip.FunctionName = 'analyze_experiment_set';
ip.addRequired('experiment_list',@isstruct);
ip.addRequired('output_dir',@isdir);
ip.addParamValue('params',struct(),@isstruct);
ip.addParamValue('load_results',false,@islogical);
ip.parse(experiment_list, output_dir, varargin{:});

% Get default parmaters
parsed_params = default_parameters(ip.Results.params);

%% Loop through experiments
%experiment_set_data = [];
experiment_counts = cell(size(ip.Results.experiment_list,2),2);
for e=1:size(ip.Results.experiment_list,2)
    experiment = ip.Results.experiment_list(e);
    fprintf('Analyzing Experiment: %s\n', experiment.name);
    exp_output_dir = [ip.Results.output_dir filesep experiment.name];
    [s,mess,messid] = mkdir(exp_output_dir); %#ok<NASGU,ASGLU>
    [experiment.spot_data experiment.cell_maps experiment.counts] = ...
        analyze_experiment(experiment.region_files, exp_output_dir, ...
        'params', parsed_params, ...
        'load_results', ip.Results.load_results, ...
        'dye_labels', experiment.dye_labels);
    %experiment_set_data = [experiment_set_data, experiment];
	experiment_counts{e,1} =  experiment.name;
    experiment_counts{e,2} = experiment.counts;
end

%% Save spot count results
save([ip.Results.output_dir filesep 'experiment_counts.mat'], 'experiment_counts')
spotCountsFile = [ip.Results.output_dir filesep 'spot_counts.tsv'];
headers = {'Experiment','Region','Cell','Cy3','Cy3.5','Cy5','Cell Phase'};
cellwrite(spotCountsFile, headers, '\t', 'wt');
for i=1:size(experiment_counts,1)
    counts = [repmat({ip.Results.experiment_list(i).name}, size(experiment_counts{i,2},1),1),  num2cell(experiment_counts{i,2})];
    cellwrite(spotCountsFile, counts, '\t', 'at');
end
