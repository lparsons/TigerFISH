% Script to read DAPI image and segment into cells using seeded watershed
% segmentation
% Based off method described on
% http://blogs.mathworks.com/steve/2006/06/02/cell-segmentation/ and in
% Cohen et al. http://www.sciencemag.org/cgi/content/abstract/322/5907/1511

function Maps = segment_cells_buds(image_input, outputdir)

%% Configuration

% Morphological open structure
morphOpenStruct = strel('disk', 6); %ones(5,5);

% Minimum cell area (pixels)
minCellSize = 600;

% Cell border thickening (pixels)
cellBorderThicken = 5;



% I = imread(filename);
%%===========================================================================================================%
% Load Image & Gets Max Projection from the Planes in Focus  %===============================================% 
%% Load image file if not already loaded
if (ischar(image_input)) 
    image = load_image(image_input);
else
    image = image_input;
end
% Parse filename
[pathstr, name, ext, versn] = fileparts(image.filename);


Std = zeros(size(image.layers,3), 1);  
for i=1:size(image.layers,3)   
    layer = image.layers(:,:,i); 
    Std(i) = std(  layer(:) ); 
end
[val ind] = max( Std ); 

Layers = image.layers(:,:,   (max(1,ind-2) : min( size(image.layers,3), ind+2) ) );
MaxProj = max(  Layers, [], 3 );
%% Scale image values
I_sc = mat2gray( MaxProj );
Maps.MaxProj = MaxProj;


%%===========================================================================================================%
% Identify Nuclei  %=========================================================================================% 

% Find nuclei - extended maxima operator can be used to identify groups of
% pixels that are significantly higher than their immediate surrounding. 
mask_em = imextendedmax( I_sc, .14);
% Cleanup nuclei using morphological close, fill,  remove small objects
% The morphological close operation is a dilation followed by an erosion,
% using the same structuring element for both operations.
mask_em = imclose(mask_em,  strel('disk', 5) );

% Fill in holes
mask_em = imfill(mask_em, 'holes');

% Remove small objects
mask_em = bwareaopen(mask_em, 20);

[Maps.nuc Maps.nucNum]  = bwlabeln( mask_em );
Maps.nucMaxproj = cell( Maps.nucNum, 1 ); 
for i=1:Maps.nucNum
	Maps.nucMaxproj{i} = Maps.MaxProj( Maps.nuc==i ); 
end



%%===========================================================================================================%
% Identify Cells  %==========================================================================================% 

%% Adaptive contrast enhancement 
% adapthisteq implements a technique called contrast-limited adaptive histogram equalization, or CLAHE.
I_eq = adapthisteq(I_sc);
%I_eq = imadjust(I);
%figure, imshow(I_eq);


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

% With n = Inf, thickens objects by adding pixels to the exterior of
% objects until doing so would result in previously unconnected objects
% being 8-connected. This option preserves the Euler number.
bw5 = bwmorph(bw4, 'thicken', cellBorderThicken);



% complement the image so that the peaks become valleys. We do this because
% we are about to apply the watershed transform, which identifies low
% points, not high points. 
I_sc_c = imcomplement(I_sc);

% modify the image so that the background pixels and the extended maxima
% pixels are forced to be the only local minima in the image. 
I_mod = imimposemin(I_sc_c, ~bw5 | mask_em);

% Compute the watershed transform
L = watershed(I_mod);
%figure, imshow(label2rgb(L));

% Get size of regions from watershed regions (cells)
S = regionprops(L, 'Area');
% Discard background (largest)
cellMap = ismember(L, find([S.Area] < max([S.Area])));

%[Maps.cells Maps.CellNum] = bwlabel( cellMap );
[Maps.cellsPerim Maps.cells Maps.CellNum] = bwboundaries(cellMap);
if 	max( Maps.cells(:) )>Maps.CellNum, 
	Maps.cells( Maps.cells>Maps.CellNum ) = 0; 
	Area = regionprops( Maps.cells, 'Area' );
	for i=1:numel(Area.Area)
		if 	Area(i).Area > 30e3
			Maps.cells( Maps.cells==i ) = 0;
		end 
	end
end
%save( 'Maps_3', 'Maps' );



%%===========================================================================================================%
% Computes probability of a cell to be budded  %=============================================================%
	
	%Custom Code for computing stats that can be returned by the function "regionprops" in the image processing toolbox 
	%1) Elipticity 
    elp.numP = zeros( Maps.CellNum, 1 ); 
	elp.axis = zeros( Maps.CellNum, 2 ); 
    elp.MaxInds = zeros( Maps.CellNum, 1 );
    elp.MinInds = zeros( Maps.CellNum, 1 ); 
    elp.val = zeros( Maps.CellNum, 1 );
	Dist = zeros( Maps.CellNum, 2 );
	Corr = zeros( Maps.CellNum, 2 );
	for i=1:Maps.CellNum
	  xy = Maps.cellsPerim{i};
	  numP = floor( size(xy,1)/2 );
	  elp.numP(i) = numP;
      
	  PerimDist = sum(  ( xy( 1:numP, : ) - xy( numP:(2*numP-1), : ) ).^2, 2 );
       [val ind] = max( PerimDist );
       elp.MaxInds(i,1) = ind;
	  MajorAxis = sqrt(  val  );
       [val ind] = min( PerimDist );
       elp.MinInds(i,1) = ind;
	  MinorAxis = sqrt( val );
	  elp.axis(i,1:2) = [MajorAxis MinorAxis]; 
      elp.val (i,1) =   sqrt( 1 - (MinorAxis / MajorAxis)^2 );
    
		clear dist1 
		prm = Maps.cellsPerim{i};
		prm = [prm; prm; prm]; 
		maxi = elp.MaxInds(i)+size(prm,1)/3;
		
		sz = round( size(prm,1)/(3*3) );
			
		dist1(:,1)= sum(  (prm( (maxi+(1:sz)), : )  - prm( (maxi-(1:sz)),  :  ) ).^2, 2 );
		
		dist1(:,2) = sum( ( prm( (maxi+numP+(1:sz)), : ) -...
							prm( (maxi+numP-(1:sz)), : )     ).^2,  2 );    

		Dist(i,:) = mean( sqrt(dist1) );  
		Corr(i,:) = corr( (1:sz)',  dist1, 'Type', 'Spearman'  );
	end
	Maps.bud.eccen=elp.val;
    Maps.bud.Dist=Dist;
	Maps.bud.Corr=Corr;
	Maps.bud.budded=zeros(size(elp.val));

	%figure, imagesc( Maps.MaxProj );
    imagesc(  (Maps.cells>0) +(Maps.nuc>0)  ); hold on 
	set( gca, 'Xtick', [], 'Ytick', [] );
	colormap( [0 0 0;  0 0 1; 1 1 1] );

     


sel = find(  sum( Corr<0.7, 2 )>0  | ...
                ( Dist(:,1) < .5*median( Dist(:,1)  ) | Dist(:,2) < .5*median( Dist(:,2)  ) ) |... 
                  elp.val( : )  >0.8   );


[MaxX MaxY ] = size( Maps.cells ); 

for i=1:numel( sel ); 
       
    j =   sel(i);
    
    if   sum( Maps.cellsPerim{ j }(:) ==1 ) || ...
         sum( Maps.cellsPerim{ j }(:,1) >= MaxX-1  ) || ... 
         sum( Maps.cellsPerim{ j }(:,2) >= MaxY-1 ),
         continue 
    end
	Maps.bud.budded( j ) = 1; 
    
    plot( Maps.cellsPerim{ j }(:,2),...
          Maps.cellsPerim{ j }(:,1), 'r', 'linewidth', 2 ); 
end	










% Output
if 0 & nargin > 1
    [s,mess,messid] = mkdir(outputdir);
    % Cell Map
    imwrite(cellMap, [outputdir filesep name '_cellmap.tif']);
    % Cell Map Colored
    cellMapRgb = label2rgb(bwlabel(cellMap), 'jet', 'k');
    imwrite(cellMapRgb, [outputdir filesep name '_cellmap_rgb.tif']);
    % Cell Perimeter
    cell_perim = bwperim(cellMap);
    cell_overlay = imoverlay(I_sc, cell_perim, [.3 .3 1]);
    imwrite(cell_overlay, [outputdir filesep name '_cellmap_overlay.tif']);
end
%figure, imshow(I), hold on
%himage = imshow(Lrgb);
%set(himage, 'AlphaData', 0.3);
%title('Segmentation superimposed transparently on original image');
%transparentCombined = getframe(get(0,'CurrentFigure'));
% merged = immerge(I, Lrgb, 0.3);
% imwrite(transparentCombined.cdata, [pathstr filesep name '_celloverlay_out.tif']);
% close;
end
