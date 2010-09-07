function [ in_focus_layer ] = best_focus_layer( multi_layer_image )
%BEST_FOCUS_LAYER Find layer with largest difference between max intensity
%and min intensity
%   Detailed explanation goes here

num_stacks = size(multi_layer_image,3);
layer_int_diff = zeros(1,num_stacks);
for l=1:num_stacks
    layer_intensities = multi_layer_image(:,:,l);
    layer_int_diff(l) = max(layer_intensities(:)) - min(layer_intensities(:));
end
[max_diff, in_focus_layer] = max(layer_int_diff);
end
