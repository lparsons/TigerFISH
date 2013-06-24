function threshold = determine_threshold(varargin)
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------

%% Parse Arguments

ip = inputParser;
ip.FunctionName = 'determine_thresholds';
ip.addRequired('out_spot_intensities',@isnumeric);
ip.addRequired('in_spot_intensities',@isnumeric);
ip.parse(varargin{:});

sorted_out_spot_intensities = sort(ip.Results.out_spot_intensities);
threshold_index = floor(size(sorted_out_spot_intensities,1)*.9);
if threshold_index > size(sorted_out_spot_intensities,1) || threshold_index == 0
    threshold = 0;
else
    threshold = sorted_out_spot_intensities(threshold_index);
end

end