function wrapper( Path, Experiment_Numbers )


warning off
addpath bin/ % /Genomics/fafner  % /Genomics/fafner/grid/users/nslavov/

if isempty( Path ), Path = 1; end
if isnumeric( Path )
    switch Path
        case{ 1, 'nslavov' }
            Path = 'fish_img/';  t = cputime;
            Path = [Path 'new/nslavov/' ];
            nm = 'my'; 
            FirstLetter = 'N';
        case{ 2, 'sandy' }
            Path = 'fish_img/'; 
            Path = [Path 'new/' ];   
            nm = 'all'; 
            FirstLetter = 'E';
    end
end



Folder.Output = nm;                           mkdir( Folder.Output );
Folder.spOut = [nm '/out'];                  mkdir( Folder.spOut );
Folder.img =   [nm '.imgs'];                mkdir( Folder.img );
Folder.hist =  [nm '.hist'];                    mkdir( Folder.hist );
Folder.PMF =   [nm '.jdis'];                 mkdir( Folder.PMF );




main_dir_list = dir( Path );



for I = 1: size(main_dir_list,1)

    if   strcmp( main_dir_list(I).name(1),  '.' ) ||...
       ~strcmp( main_dir_list(I).name(1),  FirstLetter ) 
         continue
    end
    
	
	seps = strfind( main_dir_list(I).name, '_' ); 
	Experiment_Number = str2double( main_dir_list(I).name(3:seps-1) );
    
    if ~ismember( Experiment_Number, Experiment_Numbers  ) % sum( Experiment_Number == Experiment_Numbers  ) 
        continue
    end
    
    
    
    
    
    
    Set = main_dir_list(I).name;
	subDir_list = dir( [Path Set filesep] );
    Segments = regexp(Set, '_', 'split' );
    
    File_Num = zeros(4,1);
	for ii=1:size(subDir_list,1)
         if strcmp( subDir_list(ii).name(1), '.' )
            continue
         end
	 subDir=[subDir_list(ii).name filesep];
	 file_list = dir( [Path Set filesep subDir] );

	
    for i=1:size(file_list,1)
        
        if strcmp( file_list(i).name(1), '.' )
           continue
        end
        
        if file_list(i).bytes < 50e6, continue; end 
            
        if   ~isempty(findstr(file_list(i).name, '_CCY3.tiff'))
              File_Num(1) = File_Num(1) + 1;
              cy3_file{File_Num(1)} = [subDir file_list(i).name];
              
        end
        if  ~isempty(findstr(file_list(i).name, '_CCY3.5')) 
              File_Num(2) = File_Num(2) + 1;
              cy4_file{File_Num(2)} = [subDir file_list(i).name];
        end
		if  ~isempty(findstr(file_list(i).name, '_CCY5.tiff')) 
              File_Num(3) = File_Num(3) + 1;
			  cy5_file{File_Num(3)} = [subDir file_list(i).name];
		end
        if  ~isempty(findstr(file_list(i).name, '_CDAPI.tiff')) 
              File_Num(4) = File_Num(4) + 1;
              dapi_file{File_Num(4)} = [subDir file_list(i).name];
        end
    end
    end

    if sum( File_Num ) == 0, continue, end 
    
    fprintf( '\n\n' );
    fprintf( '%s\n', Set ); 
    fprintf( 'Number of Files: %d\n', File_Num );

    if sum(File_Num==0)>=1,
        warning( 'Missing file! I will skip the field !!!' ); 
        continue
    end
    
   for reg=1:min(File_Num)
    
    %fprintf( '\nWorking on field: %d\n', reg );	
       
    %explicit check: Do the files correspond to the same region ?
    rPos_1 = max( strfind( cy3_file{reg}, 'Position ' ) )+9;
    rPos_2 = max( strfind( cy3_file{reg}, '_' ) )-1;
    
    rNum(1) = str2double( cy3_file{reg}(rPos_1:rPos_2) );
    rNum(2) = str2double( cy4_file{reg}(rPos_1:rPos_2) );
    rNum(3) = str2double( cy5_file{reg}(rPos_1:rPos_2) );
    rNum(4) = str2double( dapi_file{reg}(rPos_1:rPos_2) );
    if sum(  rNum==rNum(1)  ) ~= 4
        warning( 'Skipping files from non corresponding regions !!!' ); 
        continue
    end
    
    fprintf( '%s\t', cy3_file{reg}  );
    fprintf( '%s\t', cy4_file{reg}  );
    fprintf( '%s\t', cy5_file{reg}  );
    fprintf( '%s\n', dapi_file{reg}  );
    
    
   end
 
   fprintf( '%s\n', Set )
   fprintf( '%s\tCCY3\n', Segments{2} )
   fprintf( '%s\tCCY3.5\n', Segments{3} )
   fprintf( '%s\tCCY5\n', Segments{4} )
   fprintf( '%s\tCDAPI\n', 'DNA' )
   
    

        
   
   
   
end   
