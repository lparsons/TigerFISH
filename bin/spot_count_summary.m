function counts = spot_count_summary(varargin)

ip = inputParser;
ip.FunctionName = 'determine_thresholds';
ip.addRequired('spot_data',@isstruct);
ip.addRequired('cell_maps',@iscell);
ip.addRequired('threshold',@isstruct);
ip.parse(varargin{:});

counts = [];

dyes = fields(ip.Results.spot_data);
regions = unique(ip.Results.spot_data.(dyes{1})(:,1));

for r=1:size(regions,1)
    region_cell_map = ip.Results.cell_maps{r};
    [cell_map_labeled,num_cells] = bwlabel(region_cell_map);
    r_counts = repmat(r,num_cells+1,1);
    r_counts = horzcat(r_counts, (0:num_cells)');
    for d=1:size(dyes,1)
        dye = dyes{d};
        if (~isempty(ip.Results.spot_data.(dye)))
            region_spot_data = ip.Results.spot_data.(dye)(ip.Results.spot_data.(dye)(:,1)==r,2:7);
            dye_count = repmat(0,num_cells+1,1);
            for c=0:num_cells
                dye_count(c+1) = sum(region_spot_data(:,5)==c & region_spot_data(:,4)>ip.Results.threshold.(dye) );
            end
            r_counts = horzcat(r_counts, dye_count);
        else
            r_counts = horzcat(r_counts, repmat(0,size(r_counts,1),1));
        end
    end
    counts = vertcat(counts, r_counts);
end

% thresh_string = '';
% sep = '';
% for t=1:size(threshold,2)
%         thresh_string = [thresh_string sep num2str(threshold(t))];
%         sep = '_';
% end

% save([output_directory filesep 'spot_counts_' thresh_string '.mat'], 'counts');
% csvwrite([output_directory filesep 'spot_counts_' thresh_string '.csv'], counts);