function num_files = combine_tiffs(input_dir,output_dir,filename_re, types)
% COMBINE_TIFFS - Combine separate tiff images into stack based on filename
% 
%   num_files = combine_tiffs(input_dir,output_dir,filename_re, types)
%
%   INPUT
%       input_dir - directory with files to combine
%       output_dir - directory to save files (will be overwritten!)
%       filename_re - regex with groups for basename, dye, and z
%           default: '(?<basename>[A-Za-z0-9_\-]+?_)(?<dye>[A-Za-z0-9\.]+)_z(?<z>[0-9]{4}).tif'
%       types - types of dyes (identified by dye group of filename_re)
%           default: {'DAPI','CY3', 'CY5'}
%
%   OUTPUT
%       num_files = number of combined files created
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------

%% Config

if nargin < 2
    error('Must specify input_dir and output_dir')
end
if ~isdir(input_dir)
    error('Input dir "%s" not found', input_dir)
end
if ~isdir(output_dir)
    error('Output dir "%s" not found', output_dir)
end

if nargin < 3
    filename_re = '(?<basename>[A-Za-z0-9_\- ]+?_)(?<dye>[A-Za-z0-9\.]+)_z(?<z>[0-9]{4}).tif';
end

if nargin < 4
    types = {'DAPI','CY3','CY5'};
end


%% Get basenames
basenames = {};
all_files = dir(input_dir);
for i = 1:size(all_files,1)
    re_names = regexp(all_files(i).name, filename_re, 'names');
    if ~isempty(re_names)
        basenames = [basenames {re_names.basename}]; %#ok<AGROW>
    end
end
basenames = unique(basenames);


%% Combine images
num_files = 0;
for b = 1:size(basenames,2)
    basename = basenames{b};
    for t = 1:size(types,2)
        filebase = sprintf('%s%s', [basename, types{t}]);
        files=dir([input_dir, filesep, sprintf('%s*', filebase)]);
        disp(' ')
        disp(['filebase = ', filebase])
        for f = 1:size(files,1)
            disp(files(f).name)
            img=imread([input_dir filesep files(f).name]);
            if f == 1
                num_files = num_files + 1;
                writemode = 'overwrite';
            else
                writemode = 'append';
            end
            imwrite(img, [output_dir, filesep, filebase, '.tif'], 'tif', 'WriteMode', writemode);
        end
    end
end
