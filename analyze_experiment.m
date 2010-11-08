function [experiment_spot_data experiment_cell_maps experiment_counts] = analyze_experiment( varargin )
% analyze_experiment - analyzes each region to determine cell boundaries and spot locations and intensities
%    - determines experiment wide intensity thresholds
%    - saves cell map, spot data
%    - saves spot intensity histograms and spot location overlay images
%
%
%   [experiment_spot_data experiment_cell_maps experiment_counts] =
%       analyze_experiment( region_file_list, output_dir, [algorithm], [load_results] )
%
%       region_file_list - list of files cy3file, cy3_5file, cy5file, dapifile
%       output_dir - directory to output experiment data, histograms, etc.
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
ip.addOptional('algorithm','3D',@ischar);
ip.addOptional('load_results',false,@islogical);
ip.parse(varargin{:});

algorithms = {'3D', '2D', '2D_local'};
if (~strcmpi(ip.Results.algorithm, algorithms))
    errormsg = ['Algorithm "' ip.Results.algorithm '" not recognized.  Must be one of: ' sprintf('%s, ',algorithms{:});];
    error(errormsg)
end

%% Calculate or load results
exp_data_file = [ip.Results.output_dir filesep 'experiment_data.mat'];
if ip.Results.load_results && exist(exp_data_file, 'file')
    % Load Results
    load(exp_data_file)
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
            ip.Results.region_file_list{p,4}, reg_output_dir);
        
        experiment_spot_data.cy3 = vertcat(experiment_spot_data.cy3, horzcat(repmat(p, size(spot_data.cy3, 1),1),spot_data.cy3));
        experiment_spot_data.cy3_5 = vertcat(experiment_spot_data.cy3_5, horzcat(repmat(p, size(spot_data.cy3_5, 1),1),spot_data.cy3_5));
        experiment_spot_data.cy5 = vertcat(experiment_spot_data.cy5, horzcat(repmat(p, size(spot_data.cy5, 1),1),spot_data.cy5));
		DNA_content = [DNA_content; cell_map.DNA_content]; %MatLab is going to call vertcat anyway so I prefrer the simpler syntax
        experiment_cell_maps{p} = cell_map;
    end
    % Estimates a CDC phase for each cell based on the DNA content inferred from DAPI staining 
       %Directory and name for the plot of DNA content
       PathFileName  = [ip.Results.output_dir filesep 'DNA_Content.pdf']; 
	   DNA_Content = DNA_Content * (1/var(DNA_Content)); %normalize variance to avoid problems from treating pdf as pmf
    [cdc.phases cdc.probs] = DNA_2_cdc_phases( DNA_Content, MaxIter, PathFileName );

    %Lance, 
    %   1) give the proper value to PathFileName
    %   2) if necessary, correct my saving of cdc
       
    save(exp_data_file, 'experiment_spot_data', 'experiment_cell_maps', 'cdc' );
end



%% Determine Thresholds

dyes = fields(experiment_spot_data);
dye_color.cy3 = [0 1 0];
dye_color.cy3_5 = [1 0 0];
dye_color.cy5 = [.8 .8 .8];

%Sets a threshold of FDR
if ~exist( 'FDR_Treshold', 'var' )
    FDR_Treshold = 0.05; 
end

for d=1:size(dyes,1)
    % Determine Thresholds - 90% outside spots threshold
    dye = dyes{d};
    if (~isempty(experiment_spot_data.(dye)))
        out_spots = experiment_spot_data.(dye)(:,7)==2;
        in_spots = experiment_spot_data.(dye)(:,7)==0;
        % TODO Implement FDR calc here (need p-vals, q-vals)
        % Estimate p-values based on spots outside cells
        % - Get median, toss outliers beyond 3 sd
        %- [ycdf, xcdf] = cdfstats(int_out)
        %- p_val = interp1(x, y, intensities_in) % default is cubic
        %- [fdr q] = mafdr(p_val)
        %- threshold q values

        %Removes very bright spots outside of the cells that are likley to be artifacts and mRNAs  
        out_spots = out_spots( out_spots < 3* midian(out_spots) ); 
        %Computes the CDF for spots outside of cells
        Num = numel(out_spots);
        OUT_CDF.x = sort( out_spots );
        OUT_CDF.y = (0:1:(Num-1)) * (1/Num);      
        %Computes p values for spots inside of cells
		prob_2be_mRNA = zeros( size(prob_2be_mRNA) );
		prob_2be_mRNA( in_spots >= OUT_CDF.x(end) ) = (1 - 1/Num);
		indDimmer = find( in_spots <  OUT_CDF.x(end) );
        prob_2be_mRNA(indDimmer) = interp1( OUT_CDF.x, OUT_CDF.y, in_spots(indDimmer) );
        pvals = 1 - prob_2be_mRNA;
        qvals = mafdr( pvals );
        [val indT] = min( abs( qvals - FDR_Treshold )  );
        threshold.(dye) = in_spots( indT ); 

        %The old method
        %threshold.(dye) = determine_threshold(experiment_spot_data.(dye)(out_spots,5), experiment_spot_data.(dye)(in_spots,5));
        
        % Histogram
        histogram.(dye) = spot_intensity_histogram(experiment_spot_data.(dye)(out_spots,5), experiment_spot_data.(dye)(in_spots,5), threshold.(dye));
        title([strrep(dye, '_', '.') ' Spot Intensity Histogram'], 'FontSize', 22);
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
            spot_overlay = plot_spot_overlay(region_spot_data, threshold.(dye), size(experiment_cell_maps{reg}), 'fgcolor', dye_color.(dye), 'N', N);
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

%% Spot Count Summary
% Summary report
experiment_counts = spot_count_summary(experiment_spot_data, experiment_cell_maps, threshold);

%% Save Results
%csvwrite([ip.Results.output_dir filesep 'spot_counts.csv'], experiment_counts)
%save([exp_output_dir filesep 'experiment.mat'], 'experiment')
