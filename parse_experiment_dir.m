function experiment_list = parse_experiment_dir( varargin )
% Generate list of experiments, regions, and input image files
%
% Writes tab delimited file with the following columns:
%      Experiment, Region, Cy3_label, Cy3_file, Cy3.5_label, Cy3.5_file,
%      Cy5_label, Cy5_file, DAPI_label, DAPI_file

%% Parse Arguments
ip = inputParser;
ip.addRequired('input_dir',@isdir);
ip.addParamValue('output_filename','experiment_list.txt',@ischar);
ip.addParamValue('filemask','*',@ischar);
ip.addParamValue('region_marker','Position', @ischar);
ip.addParamValue('regex', '^(?<experiment>.+?)[\s-_]+Position[\s]+(?<region>[\d]+).+DAPI\.tiff$', @ischar);
ip.addParamValue('dyes', {'DAPI', 'CY3', 'CY3.5', 'CY5'}, @iscellstr);
ip.parse(varargin{:});

input_dir = ip.Results.input_dir;
filemask = ip.Results.filemask;
region_marker = ip.Results.region_marker;
output_filename = [ip.Results.output_filename];
[output_file fopen_message] = fopen(output_filename, 'w');
if output_file < 0
    fprintf('Unable to open file "%s"\n', output_filename)
    error(fopen_message)
end

% Get unique experiment list
exp_re = ip.Results.regex ; %['^(?<experiment>.+?)[\s-_]+' region_marker '[\s]+(?<region>[\d]+).+DAPI\.tiff$'];
exp_name_list = [];
main_dir_list = dir( [input_dir filesep filemask] );
for i = 1: size(main_dir_list,1)
    % Find DAPI files
    renames = regexp(main_dir_list(i).name, exp_re, 'names');
    if ~isempty(renames)
        exp_name_list = [exp_name_list {renames.experiment}];
    end
end
exp_name_list = unique(exp_name_list);


% Get regions for each experiment
experiment_list = [];
for i=1:size(exp_name_list,2)
    experiment.name = exp_name_list{i};
    disp(experiment.name)
    experiment.regions = [];
    exp_file_list = dir([input_dir filesep exp_name_list{i} '*C' ip.Results.dyes{1}  '.tiff']);
    for r=1:size(exp_file_list,1)
        renames = regexp(exp_file_list(r).name, exp_re, 'names');
        if ~isempty(renames) && isfield(renames,'region') && ~strcmp(renames.region, '') && strcmp(renames.experiment, exp_name_list{i})
            experiment.regions = [experiment.regions, {renames.region}];
        end
    end
    experiment.regions = unique(experiment.regions)';
    if isempty(experiment.regions)
        experiment.regions = [experiment.regions, {'1'}];
    end
    experiment.region_files = [];
    % For each region, get files
    for r=1:size(experiment.regions,1)
        if size(experiment.regions,1) > 1
            exp_reg_file_list = dir([input_dir filesep exp_name_list{i} '*' region_marker ' ' experiment.regions{r} '*.tiff']);
        else
            exp_reg_file_list = dir([input_dir filesep exp_name_list{i} '*.tiff']);
        end
        region_files = {'','','',''};
        for f=1:size(exp_reg_file_list,1)
            renames = regexp(exp_reg_file_list(f).name, exp_re, 'names');
            if ~isempty(renames) && isfield(renames, 'experiment') && strcmp(renames.experiment, exp_name_list{i})
                if  ~isempty(findstr(exp_reg_file_list(f).name, ['_C' ip.Results.dyes{1} '.tiff']))
                    region_files{4} = [input_dir filesep exp_reg_file_list(f).name];
                elseif ~isempty(findstr(exp_reg_file_list(f).name, ['_C' ip.Results.dyes{2} '.tiff']))
                    region_files{1} = [input_dir filesep exp_reg_file_list(f).name];
                elseif  ~isempty(findstr(exp_reg_file_list(f).name, ['_C' ip.Results.dyes{3} '.tiff']))
                    region_files{2} = [input_dir filesep exp_reg_file_list(f).name];
                elseif  ~isempty(findstr(exp_reg_file_list(f).name, ['_C' ip.Results.dyes{4} '.tiff']))
                    region_files{3} = [input_dir filesep exp_reg_file_list(f).name];
                end
            end
        end
        experiment.region_files = vertcat(experiment.region_files, region_files);
    end
    
    % Print to file
    for r=1:size(experiment.regions,1)
        fprintf(output_file, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', experiment.name, experiment.regions{r}, ...
            'Cy3', experiment.region_files{r,1},...
            'Cy3.5', experiment.region_files{r,2},...
            'Cy5', experiment.region_files{r,3},...
            'DAPI', experiment.region_files{r,4});
    end
    experiment_list = [experiment_list experiment];
end
end