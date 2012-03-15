function cellMap = segment_cells(image_input, second_image, params)
% Read DAPI images and segment into cells using seeded watershed
% segmentation
%
%   CELL_MAP = segment_cells( IMAGE, SECOND_IMAGE, PARAMS )
%
%       CELL_MAP = Structure with following fields:
%
%          Best_Focus_Ind: Best focus layer (int)
%                 MaxProj: Max Projection (double)
%                 cellMap: 2D Map of cells (logical)
%              cellsPerim: Cell array with boundary pixels for each cell
%                   cells: 2D Map of cells, labeled (double)
%                 CellNum: Number of cells
%         cellMap_shrink5: Cell map after shrink 5 (for DNA content est.)
%            cellsPerim_5: Cell perimeters after shrink 5
%                 cells_5: Labeled cells after shrink 5
%               CellNum_5: Number of cells after shrink 5
%                     nuc: Map of nuclei (including simply bright pixels)
%                  nucNum: Number of nuclei
%              CytoMedian: [114x1 double]
%             DNA_content: [114x1 double]
%               Cell_Size: [114x1 double]
%          Cell_2_Exclude: [114x1 double]
%                  nucPix: {1x113 cell}
%
%
%       IMAGE = Filename or image data for DAPI image
%
%       SECOND_IMAGE = Optional second image, combined with first to
%       enhance cell boundary detection (e.g. image with high
%       autoflourescence)
%
%       PARAMS = Optional struct specifying parameters to use
%           See default_paramaters.m for parameter list and default values
%
%
% Based on method described in
% http://blogs.mathworks.com/steve/2006/06/02/cell-segmentation/
% and in
% Cohen et al. http://www.sciencemag.org/cgi/content/abstract/322/5907/1511

%% Configuration
% Get default parmaters
if nargin >= 3
    parsed_params = default_parameters(params);
else
    parsed_params = default_parameters();
end

% Number of layers above a below the best infocus layer to use
layers_around_focus = parsed_params.DAPI_layers_around_focus;

% Morphological open structure
morphOpenStruct = ones(5,5);

% Cell border thickening (pixels)
cellBorderThicken = 5;


%% Load image file if not already loaded
if (ischar(image_input)) 
    image = load_image(image_input);
else
    image = image_input;
end


%% Find best focus layer and use a few layers around that
fprintf('DAPI Searching Best focus\n');
tic
DAPI_focus_layer_method = parsed_params.DAPI_focus_layer_method;
num_stacks = length(image.info);
in_focus_layer = best_focus_layer(image.layers, DAPI_focus_layer_method);

cellMap.Best_Focus_Ind = in_focus_layer;
bottom_layer = max(1,in_focus_layer-layers_around_focus);
top_layer = min(num_stacks,in_focus_layer+layers_around_focus);
max_image = max(image.layers(:,:,bottom_layer:top_layer),[],3);
cellMap.MaxProj = max_image;
% Gets the Layers closest to the best focus layer to be used for estimating
% DNA content
Layers = image.layers(:,:,  bottom_layer:top_layer  );
toc


%% Scale image values
fprintf ('Searching for nuclei\n');
tic
max_image  = max_image *( 1/median(max_image(:)) );
I_sc = mat2gray(max_image);
%I_sc = I_sc *( 1/median(I_sc(:)) );

%% Incorporate second image to better find cell boundaries
if nargin >= 2 && ~isempty(second_image)
    if (ischar(second_image)) 
        image_2 = load_image(second_image);
    else
        image_2 = second_image;
    end
    max_image_2 = max(image_2.layers(:,:,bottom_layer:top_layer),[],3);
    I_sc_2 = mat2gray(max_image_2);
    I_sc = mat2gray( I_sc + 0.2*I_sc_2 ); %*(1/median(I_sc_2(:))) );
end

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
bw4 = bwareaopen(bw3, parsed_params.minCellSize);

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
mask_em = imextendedmax(I_sc, parsed_params.DAPI_Contrast );

% Cleanup nuclei using morphological close, fill,  remove small objects
% The morphological close operation is a dilation followed by an erosion,
% using the same structuring element for both operations.
mask_em = imclose(mask_em, ones(5,5));

% Fill in holes
mask_em = imfill(mask_em, 'holes');

% Remove small objects
mask_em = bwareaopen(mask_em, 40);
toc



% Create an overlay to view
%overlay2 = imoverlay(I_eq, bw5_perim | mask_em, [.3 1 .3]);
%figure, imshow(overlay2);

fprintf('Finding cell boundaries\n')
tic
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
cellMap.cellMap = ismember(L, find([S.Area] < max([S.Area])));

% Remove spur pixels
cellMap.cellMap = bwmorph(cellMap.cellMap, 'spur', Inf);

toc


% Gets cell perimeter
[cellMap.cellsPerim cellMap.cells cellMap.CellNum] = bwboundaries(cellMap.cellMap);

% Shrinks the cells to ensure that only cytoplams is used for estimating the median DAPI auofourescence of the cytoplasm 
cellMap.cellMap_shrink5 = bwmorph( cellMap.cellMap, 'shrink', 5 );
[cellMap.cellsPerim_5 cellMap.cells_5 cellMap.CellNum_5] = bwboundaries(cellMap.cellMap_shrink5);


% Gets Coordinates and number of nuclei
fprintf('Finding nuclear pixels\n')
tic
nucs = mask_em | cellMap.MaxProj > median( cellMap.MaxProj(cellMap.cells > 0) ) * 1.5;
[cellMap.nuc cellMap.nucNum]  = bwlabeln( nucs ); 
toc



%Computes DNA contents based on DAPI intensity
fprintf('Computing DNA content\n')
tic

switch parsed_params.DAPI_Dimensions	
	case{ 2, '2D' } 
		cellMap_Layers_Cells = cellMap.cells; 
		cellMap_Layers_Cells_5 = cellMap.cells_5;
		cellMap_Layers_Nucs =  cellMap.nuc; %
		Layers = cellMap.MaxProj;	
	case{ 3, '3D' } 
		cellMap_Layers_Cells = 	  repmat( cellMap.cells,   [1 1 size(Layers,3)] ); 
		cellMap_Layers_Cells_5 =  repmat( cellMap.cells_5, [1 1 size(Layers,3)] );
		cellMap_Layers_Nucs = 	  repmat( cellMap.nuc,     [1 1 size(Layers,3)] );
    otherwise
        error('params.DAPI_Dimensions must be either 2 or 3')
end

cellMap.CytoMedian = zeros( cellMap.CellNum, 1 );
cellMap.DNA_content = zeros( cellMap.CellNum, 1 );
cellMap.Cell_Size = zeros( cellMap.CellNum, 1 );
cellMap.Cell_2_Exclude = zeros( cellMap.CellNum, 1 );

for Cell_Num=1:cellMap.CellNum
	
	% Estimates the size of the cell using the 2D projection 
	cellMap.Cell_Size(Cell_Num,1) = sum( cellMap.cells(:) == Cell_Num );
	
	% Gets the coressponding cell number in the shranken map 
	Cell_Num_5 = cellMap.cells_5( cellMap.cells==Cell_Num );
	Cell_Num_5  = Cell_Num_5( Cell_Num_5 ~= 0 ); 
	Cell_Num_5 = mode( Cell_Num_5(:) );
    
    Cell_Nuc_Pix = cellMap.nuc( cellMap.cells_5==Cell_Num_5 );
    Cell_Nuc_Pix = Cell_Nuc_Pix(Cell_Nuc_Pix>0);
    if numel( Cell_Nuc_Pix ) < 10 
       %cellMap.cells( cellMap.cells==Cell_Num ) = 0;
       cellMap.DNA_content(Cell_Num,1) = 0;
	   cellMap.Cell_2_Exclude(Cell_Num,1) = 1;
       fprintf('The %d^th Cell does not have a nucleus!!! \n', Cell_Num )
       continue
    else
        Nuc_Num = unique( Cell_Nuc_Pix ); 
    end
    
   % Gets Pixels of the nucleus that will be used for estimating DNA content
    if numel(Nuc_Num) == 1
      cellMap.nucPix{Cell_Num} = Layers( cellMap_Layers_Nucs == Nuc_Num );
    elseif numel(Nuc_Num) >= 2
        nuc1 = Layers( cellMap_Layers_Nucs == Nuc_Num(1) );
        nuc2 = Layers( cellMap_Layers_Nucs == Nuc_Num(2) );
       cellMap.nucPix{Cell_Num} = [nuc1(:); nuc2(:)];
    elseif numel(Nuc_Num) >= 3
        fprintf( 'Cell with multiple nuclei !!!\n' );
    end
    % Gets Pixels of the Cytoplasm (wirthout nucleus) that will be used for estimating DNA content      
    %Cytoplasm = Layers( cellMap_Layers_Cells_5 == Cell_Num_5 & cellMap_Layers_Nucs == 0 ); 
	Cytoplasm = Layers( cellMap_Layers_Cells == Cell_Num & cellMap_Layers_Nucs == 0 );
	
    if numel( Cytoplasm ) < 1e2
        %Gets the cell number in the map before shrinkage  
        %Cell_Num_All = mode( cellMap.cells( cellMap.cells_5==Cell_Num ) );
        %Gets the cytoplasmic pixels from the map before shrinkage 
        Cytoplasm = Layers( cellMap_Layers_Cells == Cell_Num & cellMap_Layers_Nucs == 0 );
    end
	
	if numel( Cytoplasm ) < 50 ||  numel(Cytoplasm)/numel(cellMap.nucPix{Cell_Num}) < 1.3
		cellMap.Cell_2_Exclude(Cell_Num,1) = 2; % Very large ratio of nuclear to cytoplasmic volume
	end
    %Cytoplasm = cellMap.MaxProj(  cellMap.cells == Cell_Num & cellMap.nuc == 0 );
    %Cell = cellMap.MaxProj(  cellMap.cells == Cell_Num  );
%     c1 = cellMap.cells == Cell_Num;
%     c2 = bwmorph(c1, 'thicken', 1);
%     cell_border_pixels = setdiff(find(c2), find(c1));
%     local_background = median(cellMap.MaxProj(cell_border_pixels));
    
    
	cellMap.CytoMedian(Cell_Num,1) = median(  Cytoplasm(:) );  %cellMap.nucPix
    cellMap.DNA_content(Cell_Num,1) = sum( cellMap.nucPix{Cell_Num}(:) - cellMap.CytoMedian(Cell_Num) );
end
toc




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
