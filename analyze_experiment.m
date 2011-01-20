function [experiment_spot_data experiment_cell_maps experiment_counts] = analyze_experiment( region_file_list, output_dir, varargin )
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
%   Properties
%
%   algorithm is an optional parameter that determines method of intensity measurement
%       Must be one of '3D', '2D', or '2D_local'
%       3D - Non-parametric 3D spot intensity measurement
%       2D - Uses 2D Gaussian mask with global background per image
%       2D_local - 2D Gaussian mask with local background around spot
%
%   load_results is optional parameter, if true load previous cell map and
%       spot intensity data (if it exists).  Off be default

addpath bin/

%% Parse Arguments

ip = inputParser;
ip.FunctionName = 'analyze_experiment';
ip.addRequired('region_file_list',@iscell);
ip.addRequired('output_dir',@isdir);
ip.addOptional('dye_labels',{'gene1', 'gene2', 'gene3', 'DNA'},@iscell);
ip.addParamValue('algorithm','3D',@ischar);
ip.addParamValue('load_results',false,@islogical);
ip.parse(region_file_list, output_dir, varargin{:});

algorithms = {'3D', '2D', '2D_local'};
if (~strcmpi(ip.Results.algorithm, algorithms))
    errormsg = ['Algorithm "' ip.Results.algorithm '" not recognized.  Must be one of: ' sprintf('%s, ',algorithms{:});];
    error(errormsg)
end

%% Calculate or load results
exp_data_file = [ip.Results.output_dir filesep 'experiment_data.mat'];
if ip.Results.load_results && exist(exp_data_file, 'file')
    % Load Results
    tic
    fprintf('Loading previous data...');
    load(exp_data_file)
    toc
else
    % Analyze Regions
    experiment_spot_data.cy3 = [];
    experiment_spot_data.cy3_5 = [];
    experiment_spot_data.cy5 = [];
    DNA_content = [];
    % Loop through regions
    for p=1:size(ip.Results.region_file_list,1)
        fprintf('Position number %s\n', num2str(p))
        reg_output_dir = [ip.Results.output_dir filesep 'region_' num2str(p)];
        [s,mess,messid] = mkdir(reg_output_dir);
        
        [cell_map spot_data] = analyze_region(ip.Results.region_file_list{p,1}, ...
            ip.Results.region_file_list{p,2}, ip.Results.region_file_list{p,3}, ...
            ip.Results.region_file_list{p,4}, reg_output_dir, 'debug', ip.Results.load_results);
        
        experiment_spot_data.cy3 = vertcat(experiment_spot_data.cy3, horzcat(repmat(p, size(spot_data.cy3, 1),1),spot_data.cy3));
        experiment_spot_data.cy3_5 = vertcat(experiment_spot_data.cy3_5, horzcat(repmat(p, size(spot_data.cy3_5, 1),1),spot_data.cy3_5));
        experiment_spot_data.cy5 = vertcat(experiment_spot_data.cy5, horzcat(repmat(p, size(spot_data.cy5, 1),1),spot_data.cy5));
        
        fprintf(  '%1.2f\t', cell_map.DNA_content);
        DNA_content = [DNA_content; cell_map.DNA_content]; %MatLab is going to call vertcat anyway so I prefrer the simpler syntax
        experiment_cell_maps{p} = cell_map;
    end
    % Estimates a CDC phase for each cell based on the DNA content inferred from DAPI staining
    %Directory and name for the plot of DNA content
    PathFileName  = [ip.Results.output_dir filesep 'DNA_content.pdf'];
    [cdc.phases cdc.probs] = DNA_2_cdc_phases( DNA_content, [], PathFileName );
    
    %Lance,
    %   1) give the proper value to PathFileName
    %   2) if necessary, correct my saving of cdc
    
    save(exp_data_file, 'experiment_spot_data', 'experiment_cell_maps', 'cdc' );
    % Experiment_spot_data.(dye) is matrix with 7 columns:
    %   Region, X, Y, Z, Intensity, Cell, Cell_Type
    %
    %       Cell_Type is 0 - Normal Cell (to be counted)
    %                    1 - Border Cell (ignored in counting)
    %                    2 - Background (not a cell)
end



%% Determine Thresholds

dyes = fields(experiment_spot_data);
dye_color.cy3 = [0 1 0];
dye_color.cy3_5 = [1 0 0];
dye_color.cy5 = [.8 .8 .8];

%Sets a threshold of FDR
if  ~exist( 'FDR_Treshold', 'var' )
    FDR_Treshold = 0.01;
end

for d=1:size(dyes,1)
    % Determine Thresholds - 90% outside spots threshold
    dye = dyes{d};
    if (~isempty(experiment_spot_data.(dye)))
        out_spots_ind = experiment_spot_data.(dye)(:,7)==2;
        in_spots_ind = experiment_spot_data.(dye)(:,7)==0;
        % TODO Implement FDR calc here (need p-vals, q-vals)
        % Estimate p-values based on spots outside cells
        % - Get median, toss outliers beyond 3 sd
        %- [ycdf, xcdf] = cdfstats(int_out)
        %- p_val = interp1(x, y, intensities_in) % default is cubic
        %- [fdr q] = mafdr(p_val)
        %- threshold q values
        
        %Removes very bright spots outside of the cells that are likley to be artifacts and mRNAs
        out_spots =  experiment_spot_data.(dye)(out_spots_ind,5);
        in_spots =  experiment_spot_data.(dye)(in_spots_ind,5);
        out_spots_to_keep = out_spots < 3* median(out_spots);
        out_spots = out_spots( out_spots_to_keep );
        %out_spots_ind = out_spots_ind( out_spots_to_keep );
        %Computes the CDF for spots outside of cells
        Num = sum(out_spots_to_keep);
        [OUT_CDF.x   OUT_CDF.ind ] = sort( out_spots );
        OUT_CDF.y = (0:1:(Num-1)) * (1/Num);
        
        %Computes p values for spots inside of cells
        prob_2be_mRNA.all = zeros( size(experiment_spot_data.(dye),1), 1 );
        prob_2be_mRNA.out = zeros( size(out_spots_to_keep,1), 1 );
        prob_2be_mRNA.out(out_spots_to_keep) =  OUT_CDF.y(  OUT_CDF.ind );
        prob_2be_mRNA.out(out_spots_to_keep==0) = (1 - 1/Num );
        
        prob_2be_mRNA.in = zeros( size(in_spots,1), 1 );
        prob_2be_mRNA.in( in_spots >= OUT_CDF.x(end) ) = (1 - 1/Num);
        indDimmer = (in_spots >  OUT_CDF.x(1)) & (in_spots <  OUT_CDF.x(end));
        %         [yCDF,xCDF] = cdfcalc( out_spots );
        %         OUT_CDF.x = xCDF;
        %         OUT_CDF.y = yCDF(2:end);
        nondup = [diff(OUT_CDF.x); 1] > 0;
        
        prob_2be_mRNA.in(indDimmer) = interp1( OUT_CDF.x(nondup),...
            OUT_CDF.y(nondup),...
            in_spots(indDimmer) );
        prob_2be_mRNA.all( in_spots_ind ) = prob_2be_mRNA.in;
        prob_2be_mRNA.all( out_spots_ind ) = prob_2be_mRNA.out;
        pvals = 1 - prob_2be_mRNA.in;
        
        try
            [FDR qvals pi0] = mafdr( pvals ); %fprintf( '%1.2g\n', pi0 );
            sign = find(pvals<0.1);
            [val indT] = min( abs( qvals(sign) - FDR_Treshold )  );
            threshold.(dye) = in_spots( sign(indT) );
        catch
            warning( 'FDR failed' );
            [val indT] = min( abs( pvals - 0.2 )  );
            threshold.(dye)  = in_spots( indT );
            %end
            %end
        end
        %%
        %Appending spot probabilities to the spot data matrix
        experiment_spot_data.(dye) = [experiment_spot_data.(dye)  prob_2be_mRNA.all];
        
        
        %The old method
        %threshold.(dye) = determine_threshold(experiment_spot_data.(dye)(out_spots,5), experiment_spot_data.(dye)(in_spots,5));
        
        % Histogram
        histogram.(dye) = spot_intensity_histogram(experiment_spot_data.(dye)(out_spots_ind,5), experiment_spot_data.(dye)(in_spots_ind,5), threshold.(dye));
        title([strrep(dye, '_', '.') ' Spot Intensity Histogram'], 'FontSize', 22);
        set(histogram.(dye),'PaperPositionMode','auto', 'PaperSize', [10 5], 'Units', 'inches')
        set(histogram.(dye), 'Position',  [.25 .25 9.5 4.5] );
        set(histogram.(dye),'InvertHardCopy','on');
        print(histogram.(dye), '-dpdf', [ip.Results.output_dir filesep dye '_spot_intensity_histogram.pdf'], '-r0');
        close(histogram.(dye))
        
        % Spot overlay imagesregion_spot_data
        N = 1;
        regions = unique(experiment_spot_data.(dye)(:,1));
        total_cells = 0;
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
            
            % Write cell map, colored by cell cycle prediction
            cell_map_struct = experiment_cell_maps{r};
            phases = cdc.phases(total_cells+1:cell_map_struct.CellNum+total_cells);
            total_cells = total_cells + cell_map_struct.CellNum;
            cell_map_labeled = cell_map_struct.cells;
            cell_map_image = zeros(size(cell_map_labeled));
            colors = {[1,1,1], [0,0,1], [0,1,0], [1,0,0]};
            for p=0:3
                t = ismember(cell_map_labeled, find(phases==p));
                cell_map_image = imoverlay(cell_map_image, bwperim(t), colors{p+1});
            end
            cell_map_image = imoverlay(cell_map_image, bwperim(cell_map_struct.nuc), [1 1 1]);
            imwrite(cell_map_image, [reg_output_dir filesep 'cell_map_phases.png'], 'png', 'Transparency', [0, 0, 0]);
        end
    end
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
