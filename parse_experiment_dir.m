function experiment_list = parse_experiment_dir( input_dir, filemask, region_marker )
% Generate list of experiments, regions, and input image files

% Get unique experiment list
exp_re = ['^(?<experiment>.+?)[\s-_]+' region_marker '[\s]+(?<region>[\d]+).+DAPI\.tiff$'];
exp_name_list = [];
main_dir_list = dir( [input_dir filesep filemask] );
for i = 1: size(main_dir_list,1)
    % Find DAPI files
    renames = regexp(main_dir_list(i).name, exp_re, 'names');
    if isempty(renames)
       % fprintf('Skipping %s\n', main_dir_list(i).name);
    else
        exp_name_list = [exp_name_list {renames.experiment}];
    end
end
exp_name_list = unique(exp_name_list, 'rows');


% Get regions for each experiment
experiment_list = [];
for i=1:size(exp_name_list,2)
    experiment.name = exp_name_list{i};
    experiment.regions = [];
    exp_file_list = dir([input_dir filesep exp_name_list{i} '*DAPI.tiff']);
    for r=1:size(exp_file_list,1)
        renames = regexp(exp_file_list(r).name, exp_re, 'names');
        if ~isempty(renames)
            experiment.regions = [experiment.regions, {renames.region}];

         end
    end
    experiment.regions = unique(experiment.regions);
    experiment.region_files = [];
    % For each region, get files
    for r=1:size(experiment.regions,2)
        %fprintf('Exp: %s\tRegion: %s\n', experiment.name, experiment.regions{r});
        exp_reg_file_list = dir([input_dir filesep exp_name_list{i} '*' region_marker ' ' experiment.regions{r} '*.tiff']);
        region_files = [];
        for f=1:size(exp_reg_file_list,1)
             if  ~isempty(findstr(exp_reg_file_list(f).name, '_CDAPI.tiff'))
                 region_files{4} = [input_dir filesep exp_reg_file_list(f).name];
             elseif ~isempty(findstr(exp_reg_file_list(f).name, '_CCY3.tiff'))
                 region_files{1} = [input_dir filesep exp_reg_file_list(f).name];
             elseif  ~isempty(findstr(exp_reg_file_list(f).name, '_CCY3.5.tiff'))
                 region_files{2} = [input_dir filesep exp_reg_file_list(f).name];
             elseif  ~isempty(findstr(exp_reg_file_list(f).name, '_CCY5.tiff'))
                 region_files{3} = [input_dir filesep exp_reg_file_list(f).name];
             end
             % region_files = [region_files {exp_reg_file_list(f).name}];
        end
        experiment.region_files = vertcat(experiment.region_files, region_files);
    end
    
    experiment_list = [experiment_list experiment];
end
end