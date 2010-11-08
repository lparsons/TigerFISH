function experiment_set_data = main( file_list_csv, output_dir, varargin )
% wrapper function identifies images in Path with specified experiment 
%    numbers and outputs a list of data structures to use when processing.
%    Used when experiments are in separate directories and numbered.
%    For situations where all experiments are in same directory, use 
%       'parse_experiment_dir' function instead.
%
%   [EXPERIMENT_SET] = wrapper(path, experiment_numbers, filemask, output_dir)
%       Examines 'path' for experiments listed in 'experiment_numbers' and
%       categorizes them into CY3, CY3.5, CY5, and DAPI images
%
%       EXPERIMENT_SET - list experiment data structures
%           experiment.name - name of experiment
%           experiment.regions - list of experiment regions
%           experiment.region_files - list of files used for each region
%               (,1) = Cy3_file
%               (,2) = Cy3.5_file
%               (,3) = Cy5_file
%               (,4) = DAPI_file
%


i_p = inputParser;
i_p.FunctionName = 'main';
i_p.addOptional('algorithm','3D',@ischar);
i_p.addOptional('load_results',false,@islogical);
i_p.parse(varargin{:});


experiment_set = parse_experiment_list_file(file_list_csv);
experiment_set_data = analyze_experiment_set(experiment_set, output_dir, i_p.Results.algorithm, i_p.Results.load_results);