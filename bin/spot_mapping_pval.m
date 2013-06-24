function cy = spot_mapping_pval( cy, Maps, Max_Projection, Image_Path_Name, Thresh )
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------


if ~isfield( cy, 'intensity' )
    cy.intensity = 0;
    cy.cell_counts = zeros( Maps.CellNum, 1 ); 
	cy.cellsNuc = cell( Maps.CellNum, 1 );
	cy.out.intensity=[];
    return
end

out_spots = 0;
sz = max(size(cy.intensity));
%cy.km.cells = zeros(sz,1);
cy.cell_counts = zeros( Maps.CellNum, 1 ); 
cy.cells = cell( Maps.CellNum, 1 ); 
cy.cellsNuc = cell( Maps.CellNum, 1 );
cy.sz = sz;
cy.out.intensity=[];

     k=0;
     j=0;
for i=1:sz, 

    cell_ind = Maps.cells( cy.xyz(i,1), cy.xyz(i,2) ); 


	
    if cell_ind > 0 
        j=j+1;
		cy.cells{cell_ind}(  numel(cy.cells{cell_ind})+1   ) = cy.intensity(i); 
		
		nuc_ind = Maps.nuc( cy.xyz(i,1), cy.xyz(i,2) ); 
		if nuc_ind > 0 
			cy.cellsNuc{cell_ind}(  numel(cy.cellsNuc{cell_ind})+1   ) = cy.intensity(i); 
		end
		
		
		if  cy.intensity(i) > Thresh
			k=k+1;
			cy.cell_counts(cell_ind) = ...
			cy.cell_counts(cell_ind) + 1;
			cy.in.xyz(k,:) = cy.xyz(i,:);
			cy.in.intensityTresh(k) = cy.intensity(i);
		end
		cy.in.intensity(j) = cy.intensity(i);
    else
        out_spots = out_spots + 1;
        cy.out.intensity(out_spots) = cy.intensity(i);
    end
end 
cy.out.spots_num = out_spots;


if 0 & nargin > 2 && j >=1
    close all
    imagesc( Max_Projection ); colormap( bone );
    sorted = sort( Max_Projection(:) );
    set( gca, 'Clim', [ sorted(30), sorted(end-30)] ); hold on
    
    plot( cy.in.xyz(:,2), cy.in.xyz(:,1),  'ro' );
    print( '-dpng', [ Image_Path_Name  '.png' ] ); 
    
    close all
    imagesc( Maps.cells ); hold on
    colormap( bone );
    plot( cy.in.xyz(:,2), cy.in.xyz(:,1),  'ro' );
    print( '-dpng', [ Image_Path_Name  '_.png' ] );
    %print( '-dtiff', [ Image_Path_Name  '.tif' ] );
end