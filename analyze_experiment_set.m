function experiment_counts = analyze_experiment_set( experiment_list, output_dir, varargin )
global params 
% analyze_experiment_set performs FISH image analysis to determine cell boundaries
%   and determine number of signals in each cell
%
%   experiment_set_data = analyze_experiment_set( input_dir, output_dir, [algorithm], [load_results]) 
%       finds experiments in input_dir
%
%   algorithm - optional parameter that determines method of intensity measurement
%       Must be one of '3D', '2D', or '2D_local'
%       3D - Non-parametric 3D spot intensity measurement
%       2D - Uses 2D Gaussian mask with global background per image
%       2D_local - 2D Gaussian mask with local background around spot
%
%   load_results - optional parameter, if true load previous cell map and
%       spot intensity data (if it exists).  Off be default
%
%   thresholds - optional parameter that defines thresholds for spot intensity
%       Cell array with three values, for cy3, cy3.5, and cy5
%       Default is to determine these using spots inside vs. outside of
%           cells and an FDR of 0.05
%
%
%   experiment_set_data - struct with the following fields
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
p = mfilename('fullpath');
[pathstr] = fileparts(p);
addpath([pathstr filesep 'bin'])

%% Parse Arguments
algorithms = {'3D', '2D', '2D_local'};

ip = inputParser;
ip.FunctionName = 'analyze_experiment_set';
ip.addRequired('experiment_list',@isstruct);
ip.addRequired('output_dir',@isdir);
ip.addParamValue('algorithm','2D_local',@(x)any(strcmpi(x,algorithms)));
ip.addParamValue('thresholds',{},@iscell);
ip.addParamValue('load_results',false,@islogical);
ip.parse(experiment_list, output_dir, varargin{:});

%% Loop through experiments
%experiment_set_data = [];
experiment_counts = cell(size(ip.Results.experiment_list,2),2);
for e=1:size(ip.Results.experiment_list,2)
    experiment = ip.Results.experiment_list(e);
    fprintf('Analyzing Experiment: %s\n', experiment.name);
    exp_output_dir = [ip.Results.output_dir filesep experiment.name];
    [s,mess,messid] = mkdir(exp_output_dir);
    [experiment.spot_data experiment.cell_maps experiment.counts] = ...
        analyze_experiment(experiment.region_files, exp_output_dir, ...
        'algorithm', ip.Results.algorithm, ...
        'thresholds', ip.Results.thresholds, ...
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
