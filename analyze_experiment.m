function analyze_experiment( varargin )
% analyze_experiment performs FISH image analysis to determine cell boundaries
%   and determine number of signals in each cell
%
%   analyze_experiment( region_file_list, output_dir ) 
%       region_file_list is list of files cy3file, cy3_5file, cy5file, dapifile
%
addpath bin/

%% Parse Arguments

i_p = inputParser;
i_p.FunctionName = 'analyze_region';
i_p.addRequired('region_file_list',@iscellstr);
i_p.addRequired('output_dir',@ischar);
i_p.parse(varargin{:});