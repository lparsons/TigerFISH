function combined_params = default_parameters(specified_params)
% Return Default Parameter Values, overridden by those specified

%% Default Values

% Cell Segmentation
% -----------------
% Method used to determine best focus later in DAPI images
% Allowed values: variance, intensity
default_params.DAPI_focus_layer_method = 'variance';

% Number of layers above and below best focus layer to use when segmenting
% cells
default_params.DAPI_layers_around_focus = 3;

% Use 2D or 3D data when estimating DNA content of nucleus
default_params.DAPI_Dimensions = '2D';

% Contrast (H-maxima) of nuclear pixels above cells
default_params.DAPI_Contrast = 0.137;

% Minimum cell area (pixels)
default_params.minCellSize = 40;


% Spot Detection
% --------------
% Minimum ratio of inner to outer spot intesity for each dye
% [cy3 cy3_5 cy5]
default_params.Threshold_Contrast = [1.06  1.1  1.1];

% Minimum distance between neighboring spots
% (duplicates are discarded,
% keep only the brightest)
% Distance in z is downwieghted to 25% vs. x and y
default_params.spot_merge_distance = 5;


% Spot Intensity Measurement
% --------------------------
% Algorithm used to measure spot intensity
% Allowed values: 2D_local, 2D, 3D
%   2D_local - 2D Gaussian mask with local background around spot
%   2D - Uses 2D Gaussian mask with global background per image
%   3D - Non-parametric 3D spot intensity measurement
default_params.algorithm = '2D_local';


% Manual intensity thresholding
% -----------------------------
% for Cy3, Cy3.5, Cy5 
% 3D algorithm sample: { 0.22, 0.32, 0.3 };
% 2D algorithm sample: { 750, 750, 750 };
% Use NaN to automatically compute
default_params.Threshold_Intensity = { NaN, NaN, NaN }; 



% Automatic Thresholding
% ----------------------
% ???
default_params.FDR_Threshold = 0.01;

% NULL distribution calculation method (per dye, [cy3 cy3_5 cy5])
% 1 = Deterimne a threshold based on: (1) the mode and (2) assuming the distribution of single probe intensities is symmetric
% 2 = Deterimne a threshold by assuming spots outside of cells are single probes and computing the CDF for those spots
%   Note: Probabilities are always calculated using inside vs outside spots
default_params.NULL = [1 1 1];


% Misc
% ----
% Saves detailed results (can be reused later)
default_params.Save_Detailed_Results = true;


%% Override defaults with params specified
combined_params = default_params;
if nargin > 0
    specified_param_names = fieldnames(specified_params);
    for f = 1:length(specified_param_names)
        param_name = specified_param_names{f};
        combined_params.(param_name) = specified_params.(param_name);
    end
end
