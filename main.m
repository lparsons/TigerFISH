function experiment_set_data = main( experiment_list_file, output_dir, varargin )
% main function runs fish analysis on experiments listed in file_list_csv
%
%
%   [EXPERIMENT_SET_DATA] = main(experiment_list_file, output_dir, [params])
%       Runs analysis on experiments and files listed in experiment_list_file
%       Output is generated in output_dir
%
%   INPUT
%       experiment_list_file - tab delimted file with the following fields
%       Experiment, Region, Cy3_file, Cy3.5_file, Cy5_file, DAPI_file
%
%       output_dir - directory where output is stored
%           ouput_dir/
%               experiment_set_data.mat - experiment_set_data stored as
%                   .mat file
%               spot_counts.csv - tab delimted file with summary of spot
%                   counts per cell
%                   Experiment, Position, Cell, Cy3, Cy3.5, Cy5
%               EXPERIMENT_#/
%                   experiment_data.mat - detailed spot data for all
%                       regions in the experiment
%                   DYE_intensity_histogram.pdf - histogram of spot
%                       intensities for the experiment
%                   REGION_#/
%                       cell_map.png - transparent png with outline of cell
%                           locations
%                       DYE_projection.png - transparent png with adjusted
%                           image of dye suitable for vizualization
%                       DYE_spot_image.png - transparent png with locations
%                       of identified spots circled
%                           blue - spots below threshold
%                           red - spot above thrshold
%
%   OPTIONAL PARAMETERS
%       params - optional ini file containing parameter values for analysis
%           See default_parameters.ini for list of parameters and
%           documentation
%
%       load_results - optional parameter, if true load previous cell map and
%           spot intensity data (if it exists).  Off be default
%
%   OUTPUT
%       EXPERIMENT_SET_DATA - list of experiment data structures
%
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------

ip = inputParser;
ip.FunctionName = 'main';
ip.addParamValue('ini_file','',@ischar);
ip.addParamValue('load_results',false,@islogical);
ip.parse(varargin{:});

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

experiment_set = parse_experiment_list_file(experiment_list_file);
experiment_set_data = analyze_experiment_set(experiment_set, output_dir,...
    'ini_file', ip.Results.ini_file, 'load_results', ip.Results.load_results );

disp('DONE');

end

