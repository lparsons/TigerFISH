function [cell_map spot_data] = analyze_region( varargin )
% analyze_region determines cell boundaries and outputs spot locations and
%   intensities
%
%   [CELL_MAP, SPOT_DATA] = analyze_region( cy3file, cy3_5file, cy5file,
%       dapifile )
%
%       CELL_MAP = logical matrix that identifies cells
%       SPOT_DATA = data structure with fields for cy3, cy3.5, cy5
%           each field contains a matrix with spot locations and intenties
%           Columns: x, y, z, intensity, cell_number, cell_type    
%               cell_types: 0 - normal, 1 - edge, 2 - background
%   
%
%   [CELL_MAP, SPOT_DATA] = analyze_region( cy3file, cy3_5file, cy5file,
%       dapifile, output_dir )
%
%       If output_dir is specified, transparent, enchanced projections are
%           written to the directory for each dye.

addpath bin/

%% Parse Arguments

ip = inputParser;
ip.FunctionName = 'analyze_region';
ip.addRequired('cy3file',@ischar);
ip.addRequired('cy3_5file',@ischar);
ip.addRequired('cy5file',@ischar);
ip.addRequired('dapifile',@ischar);
ip.addOptional('output_dir','',@isdir);
ip.parse(varargin{:});

%[PATHSTR,NAME,EXT,VERSN] = fileparts(file);

%% Cell Segmentation
dapi_image = load_image(ip.Results.dapifile);
cell_map = segment_cells(dapi_image);
imwrite(bwperim(cell_map), [ip.Results.output_dir filesep 'cell_map.png'], 'png', 'Transparency', 0);

%% Spot Identification
spot_data = [];

% Locate spots
cy3_image = load_image(ip.Results.cy3file);
spot_locations.cy3 = find_spots(cy3_image);
cy3_5_image = load_image(ip.Results.cy3_5file);
spot_locations.cy3_5 = find_spots(cy3_5_image);
cy5_image = load_image(ip.Results.cy5file);
spot_locations.cy5 = find_spots(cy5_image);


% Quick Fix for running both algorithms on the same spots. It has to be cleaned
spot_locations.cy3 = filter_border_spots( spot_locations.cy3 );
spot_locations.cy3_5 = filter_border_spots( spot_locations.cy3_5 );
spot_locations.cy5 = filter_border_spots( spot_locations.cy5 );


% Measure spots using D Larson's 2D Gaussian Mask algorithm
spot_data.cy3 = measure_spots(spot_locations.cy3, cy3_image, 'local');
spot_data.cy3_5 = measure_spots(spot_locations.cy3_5, cy3_5_image, 'local');
spot_data.cy5 = measure_spots(spot_locations.cy5, cy5_image, 'local');


%keyboard
% Measure spots using Nikolai's 3D Non-Parametric method 
%[spot_data.cy3 sb.cy3] = measure_spots_np(spot_locations.cy3, cy3_image);
%[spot_data.cy3_5 sb.cy3_5] = measure_spots_np(spot_locations.cy3_5, cy3_5_image);
%[spot_data.cy5 sb.cy5] = measure_spots_np(spot_locations.cy5, cy5_image);



% Merge spots
% Use original (not modified) spot locations for consistency between
% algorithms
%duplicate_threshold = 5;
%spot_data.cy3 = merge_spots(replace_locations(spot_locations.cy3, spot_data.cy3), duplicate_threshold);
%spot_data.cy3_5 = merge_spots(replace_locations(spot_locations.cy3_5, spot_data.cy3_5), duplicate_threshold);
%spot_data.cy5 = merge_spots(replace_locations(spot_locations.cy5, spot_data.cy5), duplicate_threshold);


% Spot to cell mapping
if (~isempty(spot_data.cy3)) 
    spot_data.cy3(:,5:6) = map_spots_to_cells(cell_map, spot_data.cy3(:,1:3), []);
end
if (~isempty(spot_data.cy3_5))
    spot_data.cy3_5(:,5:6) = map_spots_to_cells(cell_map, spot_data.cy3_5(:,1:3), []);
end
if (~isempty(spot_data.cy5))
    spot_data.cy5(:,5:6) = map_spots_to_cells(cell_map, spot_data.cy5(:,1:3), []);
end

%% Write colorized, transparent spot enhanced images
if (~ strcmp(ip.Results.output_dir, ''))
    % Max project DAPI image
    layers_around_focus = 3;

    max_dapi_image_adj = projected_image(dapi_image.max, dapi_image.layers, layers_around_focus);
    %[PATHSTR,NAME,EXT,VERSN] = fileparts(ip.Results.dapifile);
    B = cat(3,zeros(size(max_dapi_image_adj)),zeros(size(max_dapi_image_adj)),ones(size(max_dapi_image_adj)));
    imwrite(B, [ip.Results.output_dir filesep  'dapi_projection.png'], 'png', 'Alpha', max_dapi_image_adj);

    max_cy3_image_adj = projected_image(cy3_image.max, cy3_image.layers, Inf);
    %[PATHSTR,NAME,EXT,VERSN] = fileparts(ip.Results.cy3file);
    G = cat(3,zeros(size(max_cy3_image_adj)),ones(size(max_cy3_image_adj)),zeros(size(max_cy3_image_adj)));
    imwrite(G, [ip.Results.output_dir filesep 'cy3_projection.png'], 'png', 'Alpha', max_cy3_image_adj);

    max_cy3_5_image_adj = projected_image(cy3_5_image.max, cy3_5_image.layers, Inf);
    %[PATHSTR,NAME,EXT,VERSN] = fileparts(ip.Results.cy3_5file);
    R = cat(3,ones(size(max_cy3_5_image_adj)),zeros(size(max_cy3_5_image_adj)),zeros(size(max_cy3_5_image_adj)));
    imwrite(R, [ip.Results.output_dir filesep 'cy3_5_projection.png'], 'png', 'Alpha', max_cy3_5_image_adj);

    max_cy5_image_adj = projected_image(cy5_image.max, cy5_image.layers, Inf);
    %[PATHSTR,NAME,EXT,VERSN] = fileparts(ip.Results.cy5file);
    imwrite(ones(size(max_cy5_image_adj)), [ip.Results.output_dir filesep 'cy5_projection.png'], 'png', 'Alpha', max_cy5_image_adj);
end

end

%% Subfunctions

function original_spot_locs = replace_locations(spot_locations, spot_data)
    % Replace spot locations in spot_data with original locations
    [x y z] = ind2sub( size(spot_locations), find( spot_locations ) );
    original_spot_locs = [y x z spot_data(:,4:end)];
end

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