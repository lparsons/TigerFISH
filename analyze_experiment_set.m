function experiment_set_data = analyze_experiment_set( varargin )
% analyze_experiment_set performs FISH image analysis to determine cell boundaries
%   and determine number of signals in each cell
%
%   experiment_set_data = analyze_experiment_set( input_dir, output_dir [algorithm], [load_results]) 
%       finds experiments in input_dir
%
%   algorithm is an optional parameter that determines method of intensity measurement
%       Must be one of '3D', '2D', or '2D_local'
%       3D - Non-parametric 3D spot intensity measurement
%       2D - Uses 2D Gaussian mask with global background per image
%       2D_local - 2D Gaussian mask with local background around spot
%
%   load_results is an optional parameter, if true load previous cell map and
%       spot intensity data (if it exists).  Off be default
%
%   experiment_set_data is struct with the following fields
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
addpath bin/

%% Parse Arguments

ip = inputParser;
ip.FunctionName = 'analyze_experiment_set';
ip.addRequired('experiment_list',@isstruct);
ip.addRequired('output_dir',@isdir);
ip.addOptional('algorithm','3D',@ischar);
ip.addOptional('load_results',false,@islogical);
ip.parse(varargin{:});

algorithms = {'3D', '2D', '2D_local'};
if (~strcmpi(ip.Results.algorithm, algorithms))
    errormsg = ['Algorithm "' ip.Results.algorithm '" not recognized.  Must be one of: ' sprintf('%s, ',algorithms{:});];
    error(errormsg)
end

%% Loop through experiments
experiment_set_data = [];
parfor e=1:size(ip.Results.experiment_list,2)
    experiment = ip.Results.experiment_list(e);
    fprintf('Analyzing Experiment: %s\n', experiment.name);
    exp_output_dir = [ip.Results.output_dir filesep experiment.name];
    [s,mess,messid] = mkdir(exp_output_dir);
    [experiment.spot_data experiment.cell_maps experiment.counts] = ...
        analyze_experiment(experiment.region_files, exp_output_dir, ip.Results.algorithm, ip.Results.load_results);
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
