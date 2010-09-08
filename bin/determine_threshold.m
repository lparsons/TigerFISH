function threshold = determine_threshold(varargin)

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