function [ in_focus_layer ] = best_focus_layer_var( multi_layer_image )
%BEST_FOCUS_LAYER Find layer with largest difference between max intensity


Var = zeros(size(multi_layer_image,3), 1);  
for i=1:size(multi_layer_image,3)   
    layer = multi_layer_image(:,:,i); 
    Var(i) = var(  layer(:) ); 
end
[val ind] = max( Var ); % ind = 10; 
in_focus_layer = ind;
end
