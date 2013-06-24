function [ in_focus_layer ] = best_focus_layer( multi_layer_image, varargin)
%BEST_FOCUS_LAYER Find layer using one of two methods
% variance - largest variance in intensity (default)
% intensity - largest difference between max and min intensity
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------


methods = {'intensity', 'variance'};

ip = inputParser;
ip.FunctionName = 'analyze_experiment';
ip.addRequired('multi_layer_image',@isnumeric);
ip.addOptional('method','variance',@(x)any(strcmpi(x,methods)));
ip.parse(multi_layer_image, varargin{:});

switch ip.Results.method
    % Find layer with largest difference between max intensity
    case 'intensity'
        num_stacks = size(multi_layer_image,3);
        layer_int_diff = zeros(1,num_stacks);
        for l=1:num_stacks
            layer_intensities = multi_layer_image(:,:,l);
            layer_int_diff(l) = max(layer_intensities(:)) - min(layer_intensities(:));
        end
        [max_diff, in_focus_layer] = max(layer_int_diff);
    case 'variance'
        Var = zeros(size(multi_layer_image,3), 1);  
        for i=1:size(multi_layer_image,3)   
            layer = multi_layer_image(:,:,i); 
            Var(i) = var(  layer(:) ); 
        end
        [val ind] = max( Var ); % ind = 10; 
        in_focus_layer = ind;

end
