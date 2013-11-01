function [cell_map_struct spot_data] = analyze_region(cy3file, cy3_5file, cy5file, dapifile, varargin )
% analyze_region determines cell boundaries and outputs spot locations and
%   intensities
%
%   [CELL_MAP, SPOT_DATA] = analyze_region( cy3file, cy3_5file, cy5file, dapifile, [output_dir], 'ParamName', ParamValue, ... )
%
%   INPUT
%       output_dir - if specified, transparent, enchanced projections are
%           written to the directory for each dye.
%
%   OPTIONAL PARAMETERS
%       params - optional struct containing parameter values for analysis
%           See default_parameters.m for list of parameters and
%           documentation
%
%       load_results - optional parameter, if true load previous cell map and
%           spot intensity data (if it exists).  Off by default.
%
%   OUTPUT
%       CELL_MAP = logical matrix that identifies cells
%       SPOT_DATA = data structure with fields for cy3, cy3.5, cy5
%           each field contains a matrix with spot locations and intenties
%           Columns: x, y, z, intensity, cell_number, cell_type
%               cell_types: 0 - normal, 1 - edge, 2 - background
%
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------

p = mfilename('fullpath');
[pathstr] = fileparts(p);
addpath([pathstr filesep 'bin'])

%% Parse Arguments
algorithms = {'3D', '2D', '2D_local'};

ip = inputParser;
ip.FunctionName = mfilename;
ip.addRequired('cy3file',@(x)ischar(x) || isempty(x));
ip.addRequired('cy3_5file',@(x)ischar(x) || isempty(x));
ip.addRequired('cy5file',@(x)ischar(x) || isempty(x));
ip.addRequired('dapifile',@(x)ischar(x) && exist(x, 'file'));
ip.addOptional('output_dir','',@isdir);
ip.addParamValue('debug',false,@islogical);
ip.addParamValue('ini_file','',@ischar);
ip.parse(cy3file, cy3_5file, cy5file, dapifile, varargin{:});

% Get default parmaters
parsed_params = parse_ini(ip.Results.ini_file);


reg_data_file = [ip.Results.output_dir filesep 'region_data.mat'];
if ip.Results.debug && exist(reg_data_file, 'file')
    % Load Results
    load(reg_data_file)
else
    
    %% Cell Segmentation
    fprintf('Loading DAPI Image\n');
    tic
    dapi_image = load_image(ip.Results.dapifile);
    toc
    
    fprintf('Segmenting Cells\n');
    cell_map_struct = segment_cells(dapi_image, [], ip.Results.ini_file);
    cell_map = cell_map_struct.cellMap;
    imwrite(bwperim(cell_map), [ip.Results.output_dir filesep 'cell_map.png'], 'png', 'Transparency', 0);
    
    %% Spot Identification
    spot_data = [];
    dyes = {'cy3', 'cy3_5', 'cy5'};
    for d = 1:size(dyes,2)
        dye = dyes{d};
    contrast_thresholds = parsed_params.GetValues('spot_detection', 'Threshold_Contrast');
        
        % Skip if no image
        if (isempty(ip.Results.([dye 'file'])))
            image.(dye) = [];
            spot_locations.(dye) = [];
            spot_data.(dye) = [];
        else
            
            % Locate spots
            fprintf(['Loading ' dye ' Image\n']);
            tic
            image.(dye) = load_image(ip.Results.([dye 'file']));
            toc
            fprintf(['Locating ' dye ' Spots\n']);
            tic
            spot_locations.(dye) = find_spots(image.(dye), contrast_thresholds(d));
            toc
            
            
            % Quick Fix for running both algorithms on the same spots. It has to be cleaned/deleted
            % spot_locations.(dye) = filter_border_spots( spot_locations.(dye) );
            
            % Measure spot intensities 
			algorithm = parsed_params.GetValues('spot_measurement', 'algorithm');
            if (strcmpi(algorithm, '3D'))
                % Measure spots using Nikolai's 3D Non-Parametric method
                fprintf(['Using 3D algorithm to determine spot intensity for ' dye '\n']);
                tic
                [spot_data.(dye) sb.(dye)] = measure_spots_np(spot_locations.(dye), image.(dye));
                toc
            elseif (strcmpi(algorithm, '2D'))
                % Measure spots using D Larson's 2D Gaussian Mask
                % algorithm, with global (image) background
                fprintf(['Using 2D algorithm with global image background to determine spot intensity for ' dye '\n']);
                tic
                spot_data.(dye) = measure_spots(spot_locations.(dye), image.(dye), 'global');
                toc
            elseif (strcmpi(algorithm, '2D_local'))
                % Measure spots using D Larson's 2D Gaussian Mask
                % algorithm, with local background
                fprintf(['Using 2D algorithm with local background to determine spot intensity for ' dye '\n']);
                tic
                spot_data.(dye) = measure_spots(spot_locations.(dye), image.(dye), 'local');
                toc
            else
                errormsg = ['Algorithm "' algorithm '" not recognized.  Must be one of: ' sprintf('%s, ',algorithms{:});];
                error(errormsg)
            end
            
            % Merge nearby spots
            spot_data.(dye) = merge_spots(spot_data.(dye), parsed_params.GetValues('spot_detection', 'spot_merge_distance'));
        end
        
        % Spot to cell mapping
        fprintf(['Mapping ' dye ' spots to cells\n']);
        tic
        if (~isempty(spot_data.(dye)))
            spot_data.(dye)(:,5:6) = map_spots_to_cells(cell_map, spot_data.(dye)(:,1:3), []);
        end
        toc
    end
    
    
    %% Write colorized, transparent spot enhanced images
    fprintf('Writing colorized images\n');
    tic
    if (~ strcmp(ip.Results.output_dir, ''))
        % Max project DAPI image
        layers_around_focus = 3;
        
        max_dapi_image_adj = projected_image(dapi_image.max, dapi_image.layers, layers_around_focus);
        B = cat(3,zeros(size(max_dapi_image_adj)),zeros(size(max_dapi_image_adj)),ones(size(max_dapi_image_adj)));
        imwrite(B, [ip.Results.output_dir filesep  'dapi_projection.png'], 'png', 'Alpha', max_dapi_image_adj);
        
        if (~isempty(image.cy3)); max_cy3_image_adj = projected_image(image.cy3.max, image.cy3.layers, Inf);
        else max_cy3_image_adj = zeros(size(max_dapi_image_adj));
        end
        G = cat(3,zeros(size(max_cy3_image_adj)),ones(size(max_cy3_image_adj)),zeros(size(max_cy3_image_adj)));
        imwrite(G, [ip.Results.output_dir filesep 'cy3_projection.png'], 'png', 'Alpha', max_cy3_image_adj);
        
        if (~isempty(image.cy3_5)); max_cy3_5_image_adj = projected_image(image.cy3_5.max, image.cy3_5.layers, Inf);
        else max_cy3_5_image_adj = zeros(size(max_dapi_image_adj));
        end
        R = cat(3,ones(size(max_cy3_5_image_adj)),zeros(size(max_cy3_5_image_adj)),zeros(size(max_cy3_5_image_adj)));
        imwrite(R, [ip.Results.output_dir filesep 'cy3_5_projection.png'], 'png', 'Alpha', max_cy3_5_image_adj);
        
        if (~isempty(image.cy5)); max_cy5_image_adj = projected_image(image.cy5.max, image.cy5.layers, Inf);
        else max_cy5_image_adj = zeros(size(max_dapi_image_adj));
        end
        imwrite(ones(size(max_cy5_image_adj)), [ip.Results.output_dir filesep 'cy5_projection.png'], 'png', 'Alpha', max_cy5_image_adj);
    end
    toc
end

if ip.Results.debug
    save(reg_data_file, 'cell_map_struct', 'spot_data' );
end

end

%% Subfunctions

function max_image_adj = projected_image(image_max_project, image_layers, num_layers)
% max_image_adj - Return projection of infocus region of image, enhanced to
%   show spots more clearly
intensity_cutoff_fraction = 0.0001;
pixel_cutoff = floor(numel(image_max_project(:)) * intensity_cutoff_fraction);
b = sort(image_max_project(:));
in_focus_layer = best_focus_layer(image_layers);
bottom_layer = max(1,in_focus_layer-num_layers);
top_layer = min(size(image_layers,3),in_focus_layer+num_layers);

max_image = max(image_layers(:,:,bottom_layer:top_layer),[],3);
max_image = mat2gray(max_image, [b(pixel_cutoff), b(end-pixel_cutoff)]);
max_image_adj = imadjust(max_image,[],[0, 1],2);
end
