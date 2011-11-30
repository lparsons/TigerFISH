function [experiment_spot_data experiment_cell_maps experiment_counts] = analyze_experiment( region_file_list, output_dir, varargin )
global params
% analyze_experiment - analyzes each region to determine cell boundaries and spot locations and intensities
%    - determines experiment wide intensity thresholds
%    - saves cell map, spot data
%    - saves spot intensity histograms and spot location overlay images
%
%
%   [experiment_spot_data experiment_cell_maps experiment_counts] =
%       analyze_experiment( region_file_list, output_dir, [dye_labels], 'ParameterName', ParameterValue )
%
%       region_file_list - list of files cy3file, cy3_5file, cy5file, dapifile
%       output_dir - directory to output experiment data, histograms, etc.'
%       dye_labels - optional cell array of labels for each dye
%           defaults to {'gene1', 'gene2', 'gene3', 'DNA'}
%
%
%   Parameters
%
%   algorithm - optional parameter that determines method of intensity measurement
%       Must be one of '3D', '2D', or '2D_local'
%       3D - Non-parametric 3D spot intensity measurement
%       2D - Uses 2D Gaussian mask with global background per image
%       2D_local - 2D Gaussian mask with local background around spot
%
%   thresholds - optional parameter that defines thresholds for spot intensity
%       Cell array with three values, for cy3, cy3.5, and cy5
%       Default is to determine these using spots inside vs. outside of
%           cells and an FDR of 0.05
%
%   load_results is optional parameter, if true load previous cell map and
%       spot intensity data (if it exists).  Off be default
p = mfilename('fullpath');
[pathstr] = fileparts(p);
addpath([pathstr filesep 'bin'])


%% Parse Arguments
algorithms = {'3D', '2D', '2D_local'};

ip = inputParser;
ip.FunctionName = 'analyze_experiment';
ip.addRequired('region_file_list',@iscell);
ip.addRequired('output_dir',@isdir);
ip.addOptional('dye_labels',{'gene1', 'gene2', 'gene3', 'DNA'},@iscell);
ip.addParamValue('algorithm','3D',@(x)any(strcmpi(x,algorithms)));
ip.addParamValue('thresholds',{},@iscell);
ip.addParamValue('histogram_max',{NaN, NaN, NaN},@iscell);
ip.addParamValue('load_results',false,@islogical);
ip.parse(region_file_list, output_dir, varargin{:});

exp_data_file = [ip.Results.output_dir filesep 'experiment_data.mat'];

%% Calculate or load results
if ip.Results.load_results && exist(exp_data_file, 'file')
    % Load Results
    tic
    fprintf('Loading previous data...');
    load(exp_data_file)
    toc
else
    % Analyze Regions
    experiment_spot_data = struct('cy3', [], 'cy3_5', [], 'cy5', []);
    experiment_cell_maps = cell(size(ip.Results.region_file_list,1));
    cdc = struct();
    DNA_content1 = [];
    DNA_content = [];
    Cell_Sizes = [];
    Cell_2_Exclude = [];
    % Loop through regions
    for p=1:size(ip.Results.region_file_list,1)
        fprintf('Position number %s\n', num2str(p))
        reg_output_dir = [ip.Results.output_dir filesep 'region_' num2str(p)];
        [s,mess,messid] = mkdir(reg_output_dir); %#ok<ASGLU,NASGU>
        
        if exist(ip.Results.region_file_list{p,4}, 'file')
            [cell_map spot_data] = analyze_region(ip.Results.region_file_list{p,1}, ...
                ip.Results.region_file_list{p,2}, ip.Results.region_file_list{p,3}, ...
                ip.Results.region_file_list{p,4}, reg_output_dir, ...
                'algorithm', ip.Results.algorithm, 'debug', ip.Results.load_results);
            experiment_cell_maps{p} = cell_map;
            experiment_spot_data.cy3 = vertcat(experiment_spot_data.cy3, horzcat(repmat(p, size(spot_data.cy3, 1),1),spot_data.cy3));
            experiment_spot_data.cy3_5 = vertcat(experiment_spot_data.cy3_5, horzcat(repmat(p, size(spot_data.cy3_5, 1),1),spot_data.cy3_5));
            experiment_spot_data.cy5 = vertcat(experiment_spot_data.cy5, horzcat(repmat(p, size(spot_data.cy5, 1),1),spot_data.cy5));
            
            DNA_content1 = [DNA_content1; cell_map.DNA_content]; %MatLab is going to call vertcat anyway so I prefrer the simpler syntax
            DNA_content = [DNA_content;   cell_map.DNA_content*(1/median(cell_map.DNA_content(:) )) ];
            Cell_Sizes =  [Cell_Sizes; 	  cell_map.Cell_Size];
            Cell_2_Exclude = [Cell_2_Exclude; 	  cell_map.Cell_2_Exclude];
            
            % Estimates a CDC phase for each cell based on the DNA content inferred from DAPI staining
            % Save pdf of the plot of DNA content
            [cdc.phases cdc.probs] = DNA_2_cdc_phases( DNA_content, [] );
            cdc.DNA_content_nor = DNA_content;
            cdc.DNA_content = DNA_content1;
            cdc.Cell_Sizes = Cell_Sizes;
            cdc.Cell_2_Exclude = Cell_2_Exclude;
        else
            fprintf('\tSkipping, DAPI file not found\n')
        end
    end
    
    if params.Save_Detailed_Results
        save(exp_data_file, 'experiment_spot_data', 'experiment_cell_maps', 'cdc' );
    end
    
    % Experiment_spot_data.(dye) is matrix with 7 columns:
    %   Region, X, Y, Z, Intensity, Cell, Cell_Type
    %
    %       Cell_Type is 0 - Normal Cell (to be counted)
    %                    1 - Border Cell (ignored in counting)
    %                    2 - Background (not a cell)
end

% CDC Phases code crashes often (when DAPI images are less than perfect)
% Save first, then run this
if isfield(cdc, 'DNA_content_nor') && isfield(cdc, 'phases')
    try
        plot_cdc_phases(cdc.DNA_content_nor, cdc.phases, [ip.Results.output_dir filesep 'DNA_content.pdf']);
    catch plotError
        warning('FISHIA:DNA_Content:plottingError', 'Error plotting CDC phases: %s\n', plotError.message)
    end
end

%% Determine Thresholds, Spot Probabilities, Plot Histograms

dyes = fields(experiment_spot_data);
dye_color.cy3 = [0 1 0];
dye_color.cy3_5 = [1 0 0];
dye_color.cy5 = [.8 .8 .8];

%Sets a threshold of FDR
if  isfield(params, 'FDR_Threshold')
    FDR_Threshold = params.FDR_Threshold;
else
    FDR_Threshold = 0.01;
end
if ~isfield( params, 'NULL' )
    params.NULL = [1 1 1];
end

for d=1:size(dyes,1)
    dye = dyes{d};
    dyelabel = ip.Results.dye_labels{d};
    
    % Set Thresholds
    if ~isempty(params.Threshold_Intensity{ d } )
        threshold.(dye) = params.Threshold_Intensity{ d };
    else
        if (isempty(ip.Results.thresholds))
            threshold.dye = 0; % Default if there are no spots
        else
            threshold.(dye) = ip.Results.thresholds{d};
        end
    end
    
    if (~isempty(experiment_spot_data.(dye)))
        % Index for spots inside vs outside of cells
        out_spots_ind = experiment_spot_data.(dye)(:,7)==2;
        in_spots_ind = experiment_spot_data.(dye)(:,7)==0;
        spot_intensities = experiment_spot_data.(dye)(:,5);
        
        % Probabilities are always calculated using inside vs outside spots
        mrna_probabilies = determine_mrna_probabilities(spot_intensities, out_spots_ind, in_spots_ind, params.NULL(d) );
        % Append spot probabilities to the spot data matrix
        experiment_spot_data.(dye) = [experiment_spot_data.(dye) mrna_probabilies];
        
        % Determine Thresholds if not set
        if isempty(params.Threshold_Intensity{ d } ) && isempty(ip.Results.thresholds)
            % Using FDR and inside vs. outside spots
            threshold.(dye) = determine_threshold(spot_intensities(in_spots_ind), mrna_probabilies(in_spots_ind), FDR_Threshold);
        end
        
        % Histogram
        histogram.(dye) = spot_intensity_histogram(spot_intensities(out_spots_ind), spot_intensities(in_spots_ind), threshold.(dye), ...
            'max', ip.Results.histogram_max{d});
        title([strrep(dyelabel, '_', '.') ' Spot Intensity Histogram'], 'FontSize', 22);
        set(histogram.(dye),'PaperPositionMode','auto', 'PaperSize', [10 5], 'Units', 'inches')
        set(histogram.(dye), 'Position',  [.25 .25 9.5 4.5] );
        print(histogram.(dye), '-dpdf', [ip.Results.output_dir filesep dye '_spot_intensity_histogram.pdf'], '-r0');
        close(histogram.(dye))
        
        % Spot overlay images
        N = 1;
        regions = unique(experiment_spot_data.(dye)(:,1));
        for r=1:size(regions,1)
            reg = regions(r);
            %if (~isempty(experiment_spot_data.(dye)))
            region_spot_data = experiment_spot_data.(dye)(experiment_spot_data.(dye)(:,1)==reg,2:7);
            spot_overlay = plot_spot_overlay(region_spot_data, threshold.(dye), size(experiment_cell_maps{reg}.cellMap), 'fgcolor', dye_color.(dye), 'N', N);
            % Print figure
            reg_output_dir = [ip.Results.output_dir filesep 'region_' num2str(reg)];
            spot_overlay_filename = [reg_output_dir filesep dye '_spot_image.png'];
            screen_DPI = get(0, 'ScreenPixelsPerInch');
            print(spot_overlay, '-dpng', sprintf('-r%d', N * screen_DPI), spot_overlay_filename);
            % Load and save image with transparency
            tmp = imread(spot_overlay_filename);
            imwrite(tmp,spot_overlay_filename, 'png','Transparency', [0,0,0]);
            close(spot_overlay);
        end
    end
end

%% Write cell map, colored by cell cycle prediction
total_cells = 0;
for r=1:size(experiment_cell_maps,2)
    reg_output_dir = [ip.Results.output_dir filesep 'region_' num2str(r)];
    
    cell_map_struct = experiment_cell_maps{r};
    Region_Cells = (1:cell_map_struct.CellNum)+total_cells;
    total_cells = total_cells + cell_map_struct.CellNum;

    phases = cdc.phases(Region_Cells);
    phases(  cdc.Cell_2_Exclude( Region_Cells )>0 )= 5;

    cell_map_labeled = cell_map_struct.cells;
    cell_map_image = zeros(size(cell_map_labeled));
    colors = {[1,1,1], [0,0,1], [0,1,0], [1,0,0], [1 0.1 0.7],  [1 0.7 0.1]};
    for p=0:5
        t = ismember(cell_map_labeled, find(phases==p));
        cell_map_image = imoverlay(cell_map_image, bwperim(t), colors{p+1});
    end
    cell_map_image = imoverlay(cell_map_image, bwperim(cell_map_struct.nuc), [1 1 1]);
    imwrite(cell_map_image, [reg_output_dir filesep 'cell_map_phases.png'], 'png', 'Transparency', [0, 0, 0]);
end

%% Spot Count Summary
% Summary report
experiment_counts = spot_count_summary(experiment_spot_data, ...
    experiment_cell_maps, threshold, cdc, ...
    'output_path', ip.Results.output_dir, ...
    'dye_labels', ip.Results.dye_labels);


%% Save Results
%csvwrite([ip.Results.output_dir filesep 'spot_counts.csv'], experiment_counts)
%save([exp_output_dir filesep 'experiment.mat'], 'experiment')

end

%% Determine Threshold
% Estimate p-values based on spots outside cells
% - Get median, toss outliers beyond 3 sd
%- [ycdf, xcdf] = cdfstats(int_out)
%- p_val = interp1(x, y, intensities_in) % default is cubic
%- [fdr q] = mafdr(p_val)
%- threshold q values params.Threshold_Intensity
function threshold = determine_threshold(in_spots_intensities, in_spots_probabilites, fdr)
pvals = 1 - in_spots_probabilites;
try
    [FDR qvals pi0] = mafdr( pvals ); %fprintf( '%1.2g\n', pi0 );
    sign = find(pvals<0.05);
    [val indT] = min( abs( qvals(sign) - fdr )  );
    threshold = in_spots_intensities( sign(indT) );
catch
    warning('FISHIA:Threshold:FDRFailed', 'Unable to compute FDR using pvalues');
    [val indT] = min( abs( pvals - 0.2 )  );
    threshold = in_spots_intensities( indT );
end

if  threshold< median(in_spots_intensities)
    vals = sort( in_spots_intensities );
    threshold = vals( floor(0.08*numel(vals)) );
end
end

%% Determine mrna probabilities based on intensity values for spots inside
% vs. outside cells
function mrna_probabilities = determine_mrna_probabilities(spot_intensities, out_spots_ind, in_spots_ind, NULL_Distribution)

%Removes very bright spots outside of the cells that are likley to
%be artifacts and mRNAs
out_spots =  spot_intensities(out_spots_ind);
in_spots =  spot_intensities(in_spots_ind);
out_spots_to_keep_ind = out_spots_ind & spot_intensities <5 * median(out_spots);
out_spots_to_keep = spot_intensities(out_spots_to_keep_ind);

mrna_probabilities = zeros( size(spot_intensities, 1), 1 );
%out_spots = out_spots( out_spots_to_keep );
%out_spots_ind(~out_spots_to_keep) = 0;
%out_spots_ind = out_spots_ind( out_spots_to_keep );


switch NULL_Distribution
    case{ 1, 'IN Spots' }
        % Deterimne a threshold based on: (1) the mode and(2) assuming the distribution of single probe intensities is symmetric
        if 0
            mx = max( 3*median(in_spots) );
            interval = mx/30;
            bins = [0:interval:mx, Inf];
            [i_n i_xout] = histc(in_spots, bins);
            [val ind] = max( i_n );
            Mode = bins( ind+1 );
        else
            Mode  = median(in_spots);
        end
        NULL = [ in_spots(in_spots<Mode); 2*Mode-in_spots(in_spots<Mode)];
        Num = numel(NULL);
        [OUT_CDF.x   OUT_CDF.ind ] = sort( NULL );
        OUT_CDF.y = (0:1:(Num-1)) * (1/Num);
        
    case{ 2, 'OUT Spots' }
        %Computes the CDF for spots outside of cells
        Num = sum(out_spots_to_keep_ind);
        [OUT_CDF.x   OUT_CDF.ind ] = sort( out_spots_to_keep );
        OUT_CDF.y = (0:1:(Num-1)) * (1/Num);
        
        prob_2be_mRNA.out = zeros( size(out_spots_to_keep,1), 1 );
        prob_2be_mRNA.out =  OUT_CDF.y(  OUT_CDF.ind );
        mrna_probabilities( out_spots_to_keep_ind ) = prob_2be_mRNA.out;
        
    otherwise
        error('Unrecognized NULL_Distribution setting')
        
end

%Computes p values for spots inside cells
prob_2be_mRNA.in = zeros( size(in_spots,1), 1 );
if (~isempty(OUT_CDF.x))
    
    prob_2be_mRNA.in( in_spots >= OUT_CDF.x(end) ) = (1 - 1/Num);
    indDimmer = (in_spots >  OUT_CDF.x(1)) & (in_spots <  OUT_CDF.x(end));
    %         [yCDF,xCDF] = cdfcalc( out_spots );
    %         OUT_CDF.x = xCDF;
    %         OUT_CDF.y = yCDF(2:end);
    nondup = [diff(OUT_CDF.x); 1] > 0;
    prob_2be_mRNA.in(indDimmer) = interp1( OUT_CDF.x(nondup),...
        OUT_CDF.y(nondup),...
        in_spots(indDimmer) );
else
    warning('FISHIA:MRNA_probabilities:noOutsideSpots', 'No spots outside cells found, p values set to 0!')
end

mrna_probabilities( in_spots_ind ) = prob_2be_mRNA.in;
mrna_probabilities( out_spots_ind ) = (1 - 1/Num );
end
