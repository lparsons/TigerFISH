function image = load_image(filename)
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------

image.filename = filename;

% Read image layers
image.info = imfinfo( filename );
image.layers = zeros(image.info(1).Height, image.info(1).Width, numel(image.info));
for i=1:1:numel(image.info);
    try
        image.layers(:,:,i) = imread( filename, i, 'Info', image.info );
    catch e
        image.layers(:,:,i) = imread( filename, i);
    end
end

% Compute 16-bit max and std projections
image.max = max(image.layers, [], 3);
image.std = std(image.layers, 0, 3);