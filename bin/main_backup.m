% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------
addpath bin/
Path = '/Genomics/fafner/grid/users/nslavov/fish_img/controls/';  t = cputime;
Output_Folder = 'controls_2';  % control_ 'dual';

main_dir_list = dir( Path );
Experiment_Number = 0;
for I = 1: size(main_dir_list,1)

	if ~strcmp( main_dir_list(I).name(1), 'E' ) 
		continue
	end
	Experiment_Number = Experiment_Number+1;
	
	Set = main_dir_list(I).name;

	%Set = 'E200_ssk22-16,95,250_ssk22-1,ssk22-368_yef3';
	%Set = 'E216_om45_sur4_yef3';
	%Set = 'E218_om45_sur4_yef3';
	dir_list = dir( [Path Set filesep] );
	%dir_list = dir_list.name( ~strcmp( dir_list.name(1), '.' ) );
	mRNA_counts = [];
	Regions = {};
	reg = 0;
for Directory = 1:size(dir_list,1)
   
    %Exper = 'E200_R1_ssk22-16,95,250_ssk22-1,ssk22-368_yef3/';
    if strcmp( dir_list(Directory).name(1), '.' )
        continue
    end
    reg = reg+1;
    Exper = dir_list(Directory).name;
    
    file_list = dir( [Path Set filesep Exper] );
    for i=1:size(file_list,1)
        
        if strcmp( file_list(i).name(1), '.' )
           continue
        end
        if  ~isempty(findstr(file_list(i).name, 'cy3_')) ||...
            ~isempty(findstr(file_list(i).name, 'Cy3_'))
              cy3_file = file_list(i).name;
        end
        if  ~isempty(findstr(file_list(i).name, 'cy3p5')) || ...
            ~isempty(findstr(file_list(i).name, 'Cy3p5'))
              cy4_file = file_list(i).name;
        end
		if  ~isempty(findstr(file_list(i).name, 'cy5')) || ...
			~isempty(findstr(file_list(i).name, 'Cy5'))
			  cy5_file = file_list(i).name;
		end
        if  ~isempty(findstr(file_list(i).name, 'dapi')) || ...
            ~isempty(findstr(file_list(i).name, 'DAPI')) 
              dapi_file = file_list(i).name;
        end
    end
    
    cellMap = segment_cells( [Path Set filesep Exper filesep dapi_file ] ); %'/Image_DAPI_001.tif'
    [cy3 cy3r img3] = locSpots_2( [Path Set filesep Exper filesep cy3_file  ] ); %'/Image_Cy3_001.tif' 
    [cy4 cy4r img4] = locSpots_2( [Path Set filesep Exper filesep cy4_file  ] );  % '/Image_Cy3p5_001.tif'
	[cy5 cy5r img5] = locSpots_2( [Path Set filesep Exper filesep cy5_file  ] );

	if max(size(cy3.contrast)) <=3 || ...
       max(size(cy4.contrast)) <=3 || ...
       max(size(cy5.contrast)) <=3, continue 
    end

    %%
    % Gets a cell map
    [cellMap cell_num] = bwlabel( cellMap );
    % initializes arrays for strorying the counts
    cell.cy3.counts = zeros( cell_num, 1 );
    cell.cy4.counts = zeros( cell_num, 1 ); 
    cell.cy5.counts = zeros( cell_num, 1 );

    %% Maps Cy3 spots to Cells
    out_spots = 0;
    sz = size(cy3.km.mRNAs_ind,1);
    cy3.km.cells = zeros(sz,1);
    for ii=1:sz, 
        i = cy3.km.mRNAs_ind(ii);

        cell_ind = cellMap( cy3.xyz(i,1), cy3.xyz(i,2) ); 

        cy3.km.cells(ii) = cell_ind;

        if cell_ind > 0
            cell.cy3.counts(cell_ind) = ...
            cell.cy3.counts(cell_ind) + 1;
        else
            out_spots = out_spots + 1;
        end
    end 
    cell.cy3.out_spots = out_spots;

    %% Maps Cy4 spots to Cells
    out_spots = 0;
    sz = size(cy4.km.mRNAs_ind,1);
    cy4.km.cells = zeros(sz,1);
    for ii=1:sz, 
        i = cy4.km.mRNAs_ind(ii);

        cell_ind = cellMap( cy4.xyz(i,1), cy4.xyz(i,2) ); 

        cy4.km.cells(ii) = cell_ind;

        if cell_ind > 0
            cell.cy4.counts(cell_ind) = ...
            cell.cy4.counts(cell_ind) + 1;
        else
            out_spots = out_spots + 1;
        end
    end 
    cell.cy4.out_spots = out_spots;
    
    %% Maps Cy5 spots to Cells
    out_spots = 0;
    sz = size(cy5.km.mRNAs_ind,1);
    cy5.km.cells = zeros(sz,1);
    for ii=1:sz, 
        i = cy5.km.mRNAs_ind(ii);

        cell_ind = cellMap( cy5.xyz(i,1), cy5.xyz(i,2) ); 

        cy5.km.cells(ii) = cell_ind;

        if cell_ind > 0
            cell.cy5.counts(cell_ind) = ...
            cell.cy5.counts(cell_ind) + 1;
        else
            out_spots = out_spots + 1;
        end
    end 
    cell.cy5.out_spots = out_spots;    
    
    img.Max.Vals

    Regions(reg,:) = { cell, cy3, cy4, cy5, cy3r, cy4r, cy5r };  
    
    %ind = find( cell.cy3.counts > 0 | cell.cy4.counts > 0 );
    
    mRNA_counts = [mRNA_counts ;...
                  [cell.cy3.counts cell.cy4.counts cell.cy5.counts] ];
end
save( [Output_Folder filesep Set], 'mRNA_counts', 'Regions' );
%save( ['counts_' Set], 'mRNA_counts' );

counts(Experiment_Number,:) = { Set, mRNA_counts };
Time(Experiment_Number) =  cputime - t;
	   
save( sprintf( 'counts_%d', Experiment_Number ), 'counts', 'Time' );

delete( sprintf( 'counts_%d.mat', Experiment_Number-1 ) );

end
%%















