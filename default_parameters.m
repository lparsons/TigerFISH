function default_params = default_parameters(params)
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



%% Override defaults with params specified
if nargin > 0
    default_param_names = fieldnames(default_params);
    for f = 1:length(default_param_names)
        param_name = default_param_names{f};
        if isfield(params, param_name)
            default_params.(param_name) = params.(param_name);
        end
    end
end
