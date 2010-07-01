function image = load_image(filename)

image.filename = filename;

% Read image layers
image.info = imfinfo( filename );
image.layers = zeros(image.info(1).Height, image.info(1).Width, max(size(image.info)) );
for i=1:1:size( image.info, 1 )
    image.layers(:,:,i) = imread( filename, i );
	%image.layers_grey(:,:,i) = mat2gray( image.layers(:,:,i) ); 
end

% Compute 16-bit max and std projections
image.max = max(image.layers, [], 3);
image.std = std(image.layers, 0, 3);