function [experiment_spot_data experiment_cell_maps experiment_counts] = analyze_experiment( varargin )
% analyze_experiment performs FISH image analysis to determine cell boundaries
%   and determine number of signals in each cell
%
%   analyze_experiment( region_file_list, output_dir )
%       region_file_list is list of files cy3file, cy3_5file, cy5file, dapifile
%
addpath bin/

%% Parse Arguments

ip = inputParser;
ip.FunctionName = 'analyze_experiment';
ip.addRequired('region_file_list',@iscell);
ip.addRequired('output_dir',@isdir);
ip.parse(varargin{:});

%% Analyze Regions
experiment_spot_data.cy3 = [];
experiment_spot_data.cy3_5 = [];
experiment_spot_data.cy5 = [];
% Loop through regions
for p=1:size(ip.Results.region_file_list,1)
    fprintf('Position number %s\n', num2str(p))
    reg_output_dir = [ip.Results.output_dir filesep 'region_' num2str(p)];
    [s,mess,messid] = mkdir(reg_output_dir);
    
    [cell_map spot_data] = analyze_region(ip.Results.region_file_list{p,1}, ...
        ip.Results.region_file_list{p,2}, ip.Results.region_file_list{p,3}, ...
        ip.Results.region_file_list{p,4}, reg_output_dir);
    
    experiment_spot_data.cy3 = vertcat(experiment_spot_data.cy3, horzcat(repmat(p, size(spot_data.cy3, 1),1),spot_data.cy3));
    experiment_spot_data.cy3_5 = vertcat(experiment_spot_data.cy3_5, horzcat(repmat(p, size(spot_data.cy3_5, 1),1),spot_data.cy3_5));
    experiment_spot_data.cy5 = vertcat(experiment_spot_data.cy5, horzcat(repmat(p, size(spot_data.cy5, 1),1),spot_data.cy5));
    experiment_cell_maps{p} = cell_map;
end
%% Determine Thresholds

dyes = fields(experiment_spot_data);
for d=1:size(dyes,1)
    % Determine Thresholds - 90% outside spots threshold
    dye = dyes{d};
    if (~isempty(experiment_spot_data.(dye)))
        out_spots = experiment_spot_data.(dye)(:,7)==2;
        in_spots = experiment_spot_data.(dye)(:,7)==0;
        threshold.(dye) = determine_threshold(experiment_spot_data.(dye)(out_spots,5), experiment_spot_data.(dye)(in_spots,5));
        
        % Histogram
        histogram.(dye) = spot_intensity_histogram(experiment_spot_data.(dye)(out_spots,5), experiment_spot_data.(dye)(in_spots,5), threshold.(dye));
        title([dye ' Spot Intensity Histogram']);
        set(histogram.(dye),'PaperPositionMode','auto', 'PaperSize', [12 6])
        set(histogram.(dye), 'Position',  [200 100 1200 550] );
        print(histogram.(dye), '-dpdf', [ip.Results.output_dir filesep dye '_spot_intensity_histogram.pdf'], '-r0');
        close(histogram.(dye))
        
        % Spot overlay images
        N = 2;
        regions = unique(experiment_spot_data.(dye)(:,1));
        for r=1:regions
            %if (~isempty(experiment_spot_data.(dye)))
            region_spot_data = experiment_spot_data.(dye)(experiment_spot_data.(dye)(:,1)==r,2:7);
            spot_overlay = plot_spot_overlay(region_spot_data, threshold.(dye), size(experiment_cell_maps{r}));
            % Print figure
            reg_output_dir = [ip.Results.output_dir filesep 'region_' num2str(r)];
            spot_overlay_filename = [reg_output_dir filesep dye '_spot_image.png'];
            screen_DPI = get(0, 'ScreenPixelsPerInch');
            print(spot_overlay, '-dpng', sprintf('-r%d', N * screen_DPI), spot_overlay_filename);
            % Load and save image with transparency
            tmp = imread(spot_overlay_filename);
            imwrite(tmp,spot_overlay_filename, 'png','Transparency', [0,0,0]);
            close(spot_overlay);
        end
    end
end

%% Spot Count Summary
% Summary report
experiment_counts = spot_count_summary(experiment_spot_data, experiment_cell_maps, threshold);
