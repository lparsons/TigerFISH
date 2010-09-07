% Script to read DAPI image and segment into cells using seeded watershed
% segmentation
% Based off method described on
% http://blogs.mathworks.com/steve/2006/06/02/cell-segmentation/ and in
% Cohen et al. http://www.sciencemag.org/cgi/content/abstract/322/5907/1511

function cellMap = segment_cells(image_input)

%% Configuration
% Number of layers above a below the best infocus layer to use
layers_around_focus = 3;

% Morphological open structure
morphOpenStruct = ones(5,5);

% Minimum cell area (pixels)
minCellSize = 40;

% Cell border thickening (pixels)
cellBorderThicken = 5;


%% Load image file if not already loaded
if (ischar(image_input)) 
    image = load_image(image_input);
else
    image = image_input;
end
    
% Parse filename
[pathstr, name, ext, versn] = fileparts(image.filename);

%% Find best focus layer and use a few layers around that
num_stacks = length(image.info);
in_focus_layer = best_focus_layer(image.layers);
bottom_layer = max(1,in_focus_layer-layers_around_focus);
top_layer = min(num_stacks,in_focus_layer+layers_around_focus);
max_image = max(image.layers(:,:,bottom_layer:top_layer),[],3);


%% Scale image values
I_sc = mat2gray(max_image);

%% Adaptive contrast enhancement 
% adapthisteq implements a technique called contrast-limited adaptive
% histogram equalization, or CLAHE.
I_eq = adapthisteq(I_sc);

%% Separate cells from background
%
% Global image threshold using Otsu's method 
%   - convert greyscale to binary image
threshold = graythresh(I_eq);
%threshold2 = threshold-(threshold*.5);
bw = im2bw(I_eq, threshold);

% Cleanup thresholded image
%   Fill in holes - pixes that cannot be reached by filling in the
%   background from the edge of the image (using 4 pixel neighborhood)
bw2 = imfill(bw,'holes');

% Morphological opening using a 5x5 block.  The morphological open
% operation is an erosion followed by a dilation, using the same
% structuring element for both operations.
bw3 = imopen(bw2, morphOpenStruct);

% Morphologically open binary image (remove small objects) < 40
%   Determine connected components (8 pixel neighborhood)
%   Compute area of each component
%   Remove those below specified value
bw4 = bwareaopen(bw3, minCellSize);

% Morphological closing (dilation followed by erosion).
bw5 = bwmorph(bw4, 'close');

% With n = Inf, thickens objects by adding pixels to the exterior of
% objects until doing so would result in previously unconnected objects
% being 8-connected. This option preserves the Euler number.
bw6 = bwmorph(bw5, 'thicken', cellBorderThicken);

% Get a binary image containing only the perimeter pixels of objects in
% the input image BW1. A pixel is part of the perimeter if it is nonzero
% and it is connected to at least one zero-valued pixel. The default
% connectivity is 4.
%bw5_perim = bwperim(bw5);

% The function IMOVERLAY creates a mask-based image overlay. It takes input
% image and a binary mask, and it produces an output image whose masked
% pixels have been replaced by a specified color.
% MatLab Central -
% http://www.mathworks.com/matlabcentral/fileexchange/10502
%overlay1 = imoverlay(I_eq, bw5_perim, [.3 1 .3]);
%figure, imshow(overlay1);


%% Segment cells using Marker-based watershed segmentation
%
% Find nuclei - extended maxima operator can be used to identify groups of
% pixels that are significantly higher than their immediate surrounding. 
mask_em = imextendedmax(I_sc, .137);

% Cleanup nuclei using morphological close, fill,  remove small objects
% The morphological close operation is a dilation followed by an erosion,
% using the same structuring element for both operations.
mask_em = imclose(mask_em, ones(5,5));

% Fill in holes
mask_em = imfill(mask_em, 'holes');

% Remove small objects
mask_em = bwareaopen(mask_em, 40);

% Create an overlay to view
%overlay2 = imoverlay(I_eq, bw5_perim | mask_em, [.3 1 .3]);
%figure, imshow(overlay2);

% complement the image so that the peaks become valleys. We do this because
% we are about to apply the watershed transform, which identifies low
% points, not high points. 
I_sc_c = imcomplement(I_sc);

% modify the image so that the background pixels and the extended maxima
% pixels are forced to be the only local minima in the image. 
I_mod = imimposemin(I_sc_c, ~bw6 | mask_em);

% Compute the watershed transform
L = watershed(I_mod);
%figure, imshow(label2rgb(L));

% Get size of regions from watershed regions (cells)
S = regionprops(L, 'Area');
% Discard background (largest)
cellMap = ismember(L, find([S.Area] < max([S.Area])));

% Remove spur pixels
cellMap = bwmorph(cellMap, 'spur', Inf);

% Output
%[s,mess,messid] = mkdir(outputdir);
% Cell Map
%imwrite(cellMap, [outputdir filesep name '_cellmap.tif']);
% Cell Map Colored
%cellMapRgb = label2rgb(bwlabel(cellMap), 'jet', 'k');
%imwrite(cellMapRgb, [outputdir filesep name '_cellmap_rgb.tif']);
% Cell Perimeter
%cell_perim = bwperim(cellMap);
%cell_overlay = imoverlay(I_sc, cell_perim, [.3 .3 1]);
%imwrite(cell_overlay, [outputdir filesep name '_cellmap_overlay.tif']);
%figure, imshow(I), hold on
%himage = imshow(Lrgb);
%set(himage, 'AlphaData', 0.3);
%title('Segmentation superimposed transparently on original image');
%transparentCombined = getframe(get(0,'CurrentFigure'));
% merged = immerge(I, Lrgb, 0.3);
% imwrite(transparentCombined.cdata, [pathstr filesep name '_celloverlay_out.tif']);
% close;
end
