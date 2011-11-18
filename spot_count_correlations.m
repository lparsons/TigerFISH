%% Calculates correlation as in Silverman et al, PNAS 4/13/2010
% TODO - Output gene names for each experiment instead of 1, 2, 3
% TODO - Combine experiments when pairs of genes are the same?

%% Import data
data_dir = '/Volumes/BotLabShare/FIDO/lparsons/output/sandy_2011-11';
spot_count_file = [data_dir filesep 'spot_counts.csv'];
spot_count_import = importdata(spot_count_file, '\t', 1);

%% Setup output file
output_filename = [data_dir filesep 'gene_correlations_per_experiment.tsv'];
output_file = fopen(output_filename, 'wt');
fprintf(output_file, 'Experiment\tGene A\tGene B\tFull Data Correlation\tFull Data Corr p-value\tHigh-Signal Correlation\tHigh-Signal Corr p-value\n');
fclose(output_file);


%% Get list of experiments
tmp = spot_count_import.textdata(2:end,:);
[tmp{end+1,:}] = deal('');
u = ~all(strcmp(tmp(1:end-1,1), tmp(2:end,1)),2);
experiments = spot_count_import.textdata(u,1);

%% Generate Correlations
output_file = fopen(output_filename, 'at');
format_string = '%s\t%s\t%s\t%.2f\t%.3f\t%.2f\t%.3f\n';
for e = 1:length(experiments)
    experiment_index = strcmp(spot_count_import.textdata(2:end,1), experiments(e));
    experiment_data = spot_count_import.data(experiment_index,:);
    pairs = combnk(1:3,2);
    for p = 1:size(pairs,1)
        cell_index = experiment_data(:,2)>0;
        gene1_counts = experiment_data(cell_index,2+pairs(p,1));
        gene2_counts = experiment_data(cell_index,2+pairs(p,2));
        c = corr( gene1_counts, gene2_counts );
        
        P=10000; % Number of permutations
        corrperm = zeros(P,1);
        for i=1:P; corrperm(i) = corr(gene1_counts, gene2_counts(randperm(numel(gene2_counts)))); end;
        onesidedp = mean(corrperm>c);
        twosidedp = mean(abs(corrperm)>c);
        
        [bootstat, bootsam] = bootstrp(1000,@corr,gene1_counts,gene2_counts);
        
        hs_index = gene1_counts ~= 0 & gene2_counts ~= 0;
        hs_gene1_counts = gene1_counts(hs_index);
        hs_gene2_counts = gene2_counts(hs_index);
        hsc = corr( hs_gene1_counts, hs_gene2_counts );
        
        hs_corrperm = zeros(P,1);
        for i=1:P; hs_corrperm(i) = corr(hs_gene1_counts, hs_gene2_counts(randperm(numel(hs_gene2_counts)))); end;
        hs_onesidedp = mean(hs_corrperm>hsc);
        hs_twosidedp = mean(abs(hs_corrperm)>hsc);
        
        fprintf(output_file, format_string, ...
            experiments{e}, num2str(pairs(p,1)), num2str(pairs(p,2)), c, twosidedp, hsc, hs_twosidedp);
    end
end

%% Correlations with Cell Cycle
% Setup output file
cell_cycle_output_filename = [data_dir filesep 'gene_correlations_with_cell_cycle.tsv'];
cell_cycle_output_file = fopen(cell_cycle_output_filename, 'wt');
fprintf(cell_cycle_output_file, 'Experiment\tGene\tMean Count G1\tSD G1\tMean Count S\tSD S\tMean Count G2\tSD G2\tMean Count Other\tSD Other\n');
fclose(cell_cycle_output_file);

% Generate Cell Cycle Correlations
cell_cycle_output_file = fopen(cell_cycle_output_filename, 'at');
format_string = '%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n';
for e = 1:length(experiments)
    experiment_data_file = [data_dir filesep experiments{e} filesep 'experiment_data.mat'];
    full_experiment_data = load(experiment_data_file);
    
    experiment_index = strcmp(spot_count_import.textdata(2:end,1), experiments(e));
    experiment_data = spot_count_import.data(experiment_index,:);
    cell_index = experiment_data(:,2)>0;
    
    % phases
    % 0-> other?
    % 1-> G1
    % 2-> S
    % 3-> G2
    % 4-> other?
    % 5-> other?
    for g = 1:3
        gene_data = experiment_data(cell_index,g+2);
        fprintf(cell_cycle_output_file, format_string, ...
            experiments{e},...
            ['gene ' num2str(g)], ...
            mean(gene_data(full_experiment_data.cdc.phases==1)),...
            std(gene_data(full_experiment_data.cdc.phases==1)),...
            mean(gene_data(full_experiment_data.cdc.phases==2)),...
            std(gene_data(full_experiment_data.cdc.phases==2)),...
            mean(gene_data(full_experiment_data.cdc.phases==3)),...
            std(gene_data(full_experiment_data.cdc.phases==3)),...
            mean(gene_data(full_experiment_data.cdc.phases==0 | full_experiment_data.cdc.phases>=4)),...
            std(gene_data(full_experiment_data.cdc.phases==0 | full_experiment_data.cdc.phases>=4))...
        );
    end    
end