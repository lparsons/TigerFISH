function counts = spot_count_summary(varargin)

ip = inputParser;
ip.FunctionName = 'determine_thresholds';
ip.addRequired('spot_data',@isstruct);
ip.addRequired('cell_maps',@iscell);
ip.addRequired('threshold',@isstruct);
ip.addRequired('cdc',@isstruct);
ip.addParamValue('dye_labels',{'gene1', 'gene2', 'gene3'},@iscell)
ip.addParamValue('output_path','.',@isdir);
ip.parse(varargin{:});

counts = [];

dyes = fields(ip.Results.spot_data);
regions = [];
for d=1:size(dyes,1)
    regions = unique([regions',ip.Results.spot_data.(dyes{d})(:,1)'])';
end

Dye_probs = {};
if ~ isempty(regions)
    for r=1:size(regions,1)
        region_cell_map = ip.Results.cell_maps{r}.cellMap;
        [cell_map_labeled,num_cells] = bwlabel(region_cell_map);
        r_counts = repmat(r,num_cells+1,1);
        r_counts = horzcat(r_counts, (0:num_cells)');
        dye_probs = cell(num_cells+1,3);
        for d=1:size(dyes,1)
            dye = dyes{d};
            if (~isempty(ip.Results.spot_data.(dye)))
                region_spot_data = ip.Results.spot_data.(dye)(ip.Results.spot_data.(dye)(:,1)==r,2:8);
                dye_count = zeros(num_cells+1,1);
                for c=0:num_cells

                    Is.in_the_cell = region_spot_data(:,5)==c;
                    Is.above_threshold = region_spot_data(:,4)>ip.Results.threshold.(dye);
                    Is.yes  = Is.in_the_cell & Is.above_threshold;

                    dye_count(c+1) = sum( Is.yes );
                    dye_probs{c+1,d} = spotProb_1D( region_spot_data(Is.yes,7)  );
                end
                r_counts = horzcat(r_counts, dye_count);
            else
                r_counts = horzcat(r_counts, zeros(size(r_counts,1),1));
            end
        end
        sz = size(Dye_probs,1);
        Dye_probs((sz+1):(sz+c),:) = dye_probs(2:end,:);
        counts = vertcat(counts, r_counts);
    end

    %% Computes 2D Distributions
    N = 50;
    k=0;
    Prob2D = [];
    Y = [];
    for d1=1:size(dyes,1)
        for d2=(d1+1):size(dyes,1), k=k+1;

            % Gets info for the Joint Distributions and the Plots
            Gene1 = ip.Results.dye_labels{d1};
            Gene2 = ip.Results.dye_labels{d2};
            if (~isempty(ip.Results.spot_data.(dyes{d1})) && ~isempty(ip.Results.spot_data.(dyes{d2})))
                Prob_FileName = [ip.Results.output_path filesep 'joint_dist_prob_' Gene1 '_' Gene2 '.pdf'];
                Thresh_FileName = [ip.Results.output_path filesep 'joint_dist_thresh_' Gene1 '_' Gene2 '.pdf'];

                % Probabilistic
                Prob2D{k} = zeros( N );
                for i=1:size(Dye_probs,1)
                    Prob = Dye_probs{i,d1}(:) * Dye_probs{i,d2}(:)';

                    [sz1 sz2] = size( Prob );
                    Prob2D{k}(1:min(sz1,N), 1:min(sz2,N)) = Prob2D{k}(1:min(sz1,N), 1:min(sz2,N))  + Prob(1:min(sz1,N), 1:min(sz2,N));
                end
                Y.prob = jointDist_probs( Prob2D{k}, Gene1, Gene2, Prob_FileName );
                joint_dist_prob = Prob2D{k};
                save([ip.Results.output_path filesep 'joint_dist_prob_' Gene1 '_' Gene2 '.mat'], 'joint_dist_prob');
                %csvwrite([ip.Results.output_path filesep 'joint_dist_thresh_'Gene1 '_' Gene2 '.csv'], Y.threshold);

                % Deterministic
                Y.threshold = jointDist( counts( :, d1+2 ), counts( :, d2+2 ), Gene1, Gene2, Thresh_FileName );
                save([ip.Results.output_path filesep 'joint_dist_thresh_' Gene1 '_' Gene2 '.mat'], 'Y');
                %csvwrite([ip.Results.output_path filesep 'joint_dist_prob_' Gene1 '_' Gene2 '.csv'], Prob2D{k}/sum(Prob2D{k}(:)));
            end

        end
    end
else % No region data
    counts = [];
    Y = [];
    Prob2D = [];
end
    cdc = ip.Results.cdc;
save( [ip.Results.output_path filesep 'endCountsProbs.mat'], 'counts', 'Y', 'Prob2D', 'cdc' );


% thresh_string = '';
% sep = '';
% for t=1:size(threshold,2)
%         thresh_string = [thresh_string sep num2str(threshold(t))];
%         sep = '_';
% end

% save([output_directory filesep 'spot_counts_' thresh_string '.mat'], 'counts');
% csvwrite([output_directory filesep 'spot_counts_' thresh_string '.csv'], counts);
