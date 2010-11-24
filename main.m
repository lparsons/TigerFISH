function experiment_set_data = main( experiment_list_file, output_dir, varargin )
% main function runs fish analysis on experiments listed in file_list_csv
%
%
%   [EXPERIMENT_SET_DATA] = main(experiment_list_file, output_dir)
%       Runs analysis on experiments and files listed in experiment_list_file
%       Output is generated in output_dir
%
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
%       
%       EXPERIMENT_SET_DATA - list of experiment data structuress
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
i_p.addParamValue('algorithm','3D',@ischar);
i_p.addParamValue('load_results',false,@islogical);
i_p.parse(varargin{:});


experiment_set = parse_experiment_list_file(experiment_list_file);
experiment_set_data = analyze_experiment_set(experiment_set, output_dir, i_p.Results.algorithm, i_p.Results.load_results);