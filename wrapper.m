function experiment_set_data = wrapper( Path, Experiment_Numbers, varargin )
% wrapper function identifies images in Path with specified experiment 
%    numbers and outputs a list of data structures to use when processing.
%    Used when experiments are in separate directories and numbered.
%    For situations where all experiments are in same directory, use 
%       'parse_experiment_dir' function instead.
%
%   [EXPERIMENT_SET] = wrapper(path, experiment_numbers, filemask, output_dir)
%       Examines 'path' for experiments listed in 'experiment_numbers' and
%       categorizes them into CY3, CY3.5, CY5, and DAPI images
%
%       EXPERIMENT_SET - list experiment data structures
%           experiment.name - name of experiment
%           experiment.regions - list of experiment regions
%           experiment.region_files - list of files used for each region
%               (,1) = Cy3_file
%               (,2) = Cy3.5_file
%               (,3) = Cy5_file
%               (,4) = DAPI_file
%

addpath bin/

%% Parse Arguments

i_p = inputParser;
i_p.FunctionName = 'wrapper';
i_p.addOptional('filemask','*',@ischar);
i_p.addOptional('output_dir','output',@ischar);
i_p.addOptional('algorithm','3D',@ischar);
i_p.addOptional('load_results',false,@islogical);
i_p.parse(varargin{:});
nm = i_p.Results.output_dir;
filemask = i_p.Results.filemask;

if isempty( Path ), Path = 1; end
if isnumeric( Path )
    switch Path
        case{ 1, 'nslavov' }
            Path = '/Genomics/grid/users/nslavov/locSpot/fish_img/';  % t = cputime;
            Path = [Path 'new/nslavov/' ];
            nm = 'my';
            filemask = 'N*';
        case{ 2, 'sandy' }
            Path = '/Genomics/grid/users/nslavov/locSpot/fish_img/';
            Path = [Path 'new/' ];
            nm = 'all';
            filemask = 'E*';
    end
end



Folder.Output = nm;          mkdir( Folder.Output );
Folder.spOut = [nm '/out'];  mkdir( Folder.spOut );
Folder.img =   [nm '.imgs']; mkdir( Folder.img );
Folder.hist =  [nm '.hist']; mkdir( Folder.hist );
Folder.PMF =   [nm '.jdis']; mkdir( Folder.PMF );


%% Open output file
output_filename = [Folder.Output filesep 'experiment_list.txt'];
output_file = fopen(output_filename, 'w');


%% Loop through files in directory
main_dir_list = dir( [Path filesep filemask] );
for I = 1: size(main_dir_list,1)
    
    if   strcmp( main_dir_list(I).name(1),  '.' )
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
    
    fprintf('\n\n' );
    fprintf('%s\n', Set );
    %fprintf( 'Number of Files: %d\n', File_Num );
    
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
        
        fprintf( output_file, '%s\t', Set);
        fprintf( output_file, '%s\t', num2str(reg) );
        fprintf( output_file, '%s\t', [Path filesep Set filesep cy3_file{reg}]  );
        fprintf( output_file, '%s\t', [Path filesep Set filesep cy4_file{reg}]  );
        fprintf( output_file, '%s\t', [Path filesep Set filesep cy5_file{reg}]  );
        fprintf( output_file, '%s\n', [Path filesep Set filesep dapi_file{reg}]  );
        
        
    end
    
    %fprintf( '%s\n', Set );
    fprintf( '%s\tCCY3\n', Segments{2} );
    fprintf( '%s\tCCY3.5\n', Segments{3} );
    fprintf( '%s\tCCY5\n', Segments{4} );
    fprintf( '%s\tCDAPI\n', 'DNA' );
    
end
fclose(output_file);

% Read output_file and parse into experiment_set datastructure
% TODO Refactor this to create data structure directly
% experiment_set = parse_experiment_list_file(output_filename);

%% Run analysis
% experiment_set_data = analyze_experiment_set(experiment_set, nm, i_p.Results.algorithm, i_p.Results.load_results);
end