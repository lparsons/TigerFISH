function generate_experiment_list( path, varargin )
% wrapper function identifies images in 'path' with specified experiment
%    numbers and outputs a tab delimited file of experiment names, dye
%    labels and file paths.
%
%    Used when experiments are in separate directories and numbered.
%    For situations where all experiments are in same directory, use
%       'parse_experiment_dir' function instead.
%
%   generate_experiment_list(path, experiment_numbers, output_filename, filemask)
%       Examines 'path' for experiments listed in 'experiment_numbers' and
%       categorizes them into CY3, CY3.5, CY5, and DAPI images
%
%   writes tab delimited file with the following columns:
%      Experiment, Region, Cy3_label, Cy3_file, Cy3.5_label, Cy3.5_file, Cy5_label, Cy5_file, DAPI_label, DAPI_file

addpath bin/

%% Parse Arguments

ip = inputParser;
ip.FunctionName = 'generate_experiment_list';
ip.addOptional('experiment_numbers',1:1000,@isnumeric);
ip.addOptional('output_filename','experiment_list.txt',@ischar);
ip.addOptional('filemask','*',@ischar);
ip.addOptional('load_results',false,@islogical);
ip.addParamValue('nmiss',0,@isnumeric);
ip.parse(varargin{:});

exp_expression = ...
    '(?<experiment_letters>[a-zA-Z]+)(?<experiment_number>\d+)_+(?<cy3>[a-zA-Z0-9]+)_+(?<cy3p5>[a-zA-Z0-9]+)_+(?<cy5>[a-zA-Z0-9]+)_*(?<condition>.*)*';

%% Open output file
output_filename = ip.Results.output_filename;
output_file = fopen(output_filename, 'w');


%% Loop through files in directory
main_dir_list = dir( [path filesep ip.Results.filemask] );
for I = 1: size(main_dir_list,1)
    if   strcmp( main_dir_list(I).name(1),  '.' )
        continue
    end
    
    experiment_number = -1;
    name = main_dir_list(I).name;
    exp_fields = regexp(name, exp_expression, 'names');
    if (~isempty(exp_fields))
        experiment_number = str2double( exp_fields.experiment_number );
    end
    if ~ismember( experiment_number, ip.Results.experiment_numbers )
        continue
    end
    
    Set = main_dir_list(I).name;
    subDir_list = dir( [path filesep Set filesep] );
    Segments = regexp(Set, '_', 'split' );
    
    File_Num = zeros(4,1);
    for ii=1:size(subDir_list,1)
        if strcmp( subDir_list(ii).name(1), '.' )
            continue
        end
        subDir=[subDir_list(ii).name filesep];
        file_list = dir( [path filesep Set filesep subDir] );
        
        
        for i=1:size(file_list,1)
            
            if strcmp( file_list(i).name(1), '.' )
                continue
            end
            
            if file_list(i).bytes < 50e6, continue; end
            
            if   ~isempty(strfind(file_list(i).name, '_CCY3.tiff'))
                File_Num(1) = File_Num(1) + 1;
                cy3_file{File_Num(1)} = [subDir file_list(i).name];
                
            end
            if  ~isempty(strfind(file_list(i).name, '_CCY3.5'))
                File_Num(2) = File_Num(2) + 1;
                cy4_file{File_Num(2)} = [subDir file_list(i).name];
            end
            if  ~isempty(strfind(file_list(i).name, '_CCY5.tiff'))
                File_Num(3) = File_Num(3) + 1;
                cy5_file{File_Num(3)} = [subDir file_list(i).name];
            end
            if  ~isempty(strfind(file_list(i).name, '_CDAPI.tiff'))
                File_Num(4) = File_Num(4) + 1;
                dapi_file{File_Num(4)} = [subDir file_list(i).name];
            end
        end
    end
    
    if sum( File_Num ) == 0, continue, end
    
    fprintf('\n%s\n', Set );
    %fprintf( 'Number of Files: %d\n', File_Num );
    
    if sum(File_Num==0) > ip.Results.nmiss,
        warning( 'FISHIA:experiment_list:missingFile', 'Missing file! I will skip the field !!!' );
        continue
    end
    
%    if isfield( params, 'Region_num' )
%        Region_num = params.Region_num;
%    else
        Region_num = mode(File_Num);
%    end
    for reg=1: Region_num
        
        %fprintf( '\nWorking on field: %d\n', reg );
        
        %explicit check: Do the files correspond to the same region ?
        rPos_1 = max( strfind( cy3_file{reg}, 'Position ' ) )+9;
        rPos_2 = max( strfind( cy3_file{reg}, '_' ) )-1;
        
        rNum(1) = str2double( cy3_file{reg}(rPos_1:rPos_2) );
        rNum(2) = str2double( cy4_file{reg}(rPos_1:rPos_2) );
        rNum(3) = str2double( cy5_file{reg}(rPos_1:rPos_2) );
        rNum(4) = str2double( dapi_file{reg}(rPos_1:rPos_2) );
        if sum(  rNum==rNum(1)  ) < (4-ip.Results.nmiss)
            warning( 'FISHIA:experiment_list:nonMatchingRegions', 'Skipping files from non corresponding regions !!!' );
            continue
        end
        
        
        
        fprintf( output_file, '%s\t', Set);
        fprintf( output_file, '%s\t', num2str(reg) );
        fprintf( output_file, '%s\t', Segments{2}  );
        fprintf( output_file, '%s\t', [path filesep Set filesep cy3_file{reg}]  );
        fprintf( output_file, '%s\t', Segments{3}  );
        fprintf( output_file, '%s\t', [path filesep Set filesep cy4_file{reg}]  );
        if  sum(  rNum==rNum(1)  ) >= (4-ip.Results.nmiss)
            fprintf( output_file, '%s\t', Segments{4}  );
            fprintf( output_file, '%s\t', [path filesep Set filesep cy5_file{reg}]  );
        end
        fprintf( output_file, '%s\t', 'DAPI'  );
        fprintf( output_file, '%s\n', [path filesep Set filesep dapi_file{reg}]  );
        
        
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