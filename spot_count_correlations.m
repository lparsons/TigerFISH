function [correlations, experiments, headers] = spot_count_correlations(spotCountFile, varargin)
%SPOT_COUNT_CORRELATIONS -  Calculates gene correlations and p-values as in Silverman et al, PNAS 13-Apr-2010
%
%Reports two-sided Pearson correlations between spot counts of gene pairs
%using all of the data (full data) and using only cells where both genes
%were detected at least once (high signal).
%
%p-values are calculated using random permutations (default of 10000) of
%the data (without replacement) and reported as the fraction of the random
%correlations whose absolute values are greater than correlation calculated
%for the actual data
%
% Syntax:  [correlations, experiments, headers] = spot_count_correlations(spotCountFile)
%
% Inputs:
%    spotCountFile - Tab separated file containing output from
%    analyze_experiment_set
%        columns  = {'Experiment','Region','Cell','Cy3','Cy3.5','Cy5'}
%
% Outputs:
%    correlations - Array containing gene/gene correlations
%        for each pair of genes for each experiment
%    experiments - Cell array of experiment names for the rows of
%        'correlaltions'
%    headers - Cell array of headers for the columns of 'correlaltions'
%
% Example:
%    [correlations, experiments, headers] = spot_count_correlations('spot_counts.tsv');
%    correlationsCellArray = vertcat(headers,horzcat(experiments,num2cell(correlations)));
%    cellwrite('correlations.tsv', correlationsCellArray, '\t', 'wt');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------

% Author: Lance Parsons
% Lewis-Sigler Institute for Integrative Genomics, Princeton University
% email: lparsons@princeton.edu
% Website: http://www.linkedin.com/in/lparsons
% Nov 2011; Last revision: 01-Dec-2011

%------------- BEGIN CODE --------------

ip = inputParser;
ip.FunctionName = mfilename;
ip.addRequired('spotCountFile',@(x)exist(x, 'file'));
ip.addParamValue('permutations',10000,@(x)x-round(x)==0);
ip.parse(spotCountFile, varargin{:});

spot_count_import = importdata(spotCountFile, '\t', 1);

headers = {'Experiment', 'Gene A', 'Gene B', 'Full Data Correlation', ...
    'Full Data Corr p-value', 'High-Signal Correlation', ...
    'High-Signal Corr p-value'};

% %% Setup output file
% output_filename = [dataDir filesep 'gene_correlations_per_experiment.tsv'];
% output_file = fopen(output_filename, 'wt');
% fprintf(output_file, 'Experiment\tGene A\tGene B\tFull Data Correlation\tFull Data Corr p-value\tHigh-Signal Correlation\tHigh-Signal Corr p-value\n');
% fclose(output_file);


%% Get list of experiments
tmp = spot_count_import.textdata(2:end,:);
[tmp{end+1,:}] = deal('');
u = ~all(strcmp(tmp(1:end-1,1), tmp(2:end,1)),2);
experiment_list = spot_count_import.textdata(u,1);

%% Generate Correlations
correlations = zeros(length(experiment_list)*3,length(headers)-1);
experiments = cell(length(experiment_list)*3,1);
for e = 1:length(experiment_list)
    experiment_index = strcmp(spot_count_import.textdata(2:end,1), experiment_list(e));
    experiment_data = spot_count_import.data(experiment_index,:);
    pairs = combnk(1:3,2);
    for p = 1:size(pairs,1)
        cell_index = experiment_data(:,2)>0;
        gene1_counts = experiment_data(cell_index,2+pairs(p,1));
        gene2_counts = experiment_data(cell_index,2+pairs(p,2));
        [fullDataCorrelation, fullDataPvalue] = ...
            computeCorrelationAndPvalue(gene1_counts, gene2_counts, ...
            ip.Results.permutations);
        
        hs_index = gene1_counts ~= 0 & gene2_counts ~= 0;
        hs_gene1_counts = gene1_counts(hs_index);
        hs_gene2_counts = gene2_counts(hs_index);
        [highSignalCorrelation, highSignalPvalue] = ...
            computeCorrelationAndPvalue(hs_gene1_counts, hs_gene2_counts, ...
            ip.Results.permutations);
        
        correlations(((e-1)*3)+p,:) = ...
            [pairs(p,1), pairs(p,2), ...
            fullDataCorrelation, fullDataPvalue, ...
            highSignalCorrelation, highSignalPvalue];
        experiments{((e-1)*3)+p} = experiment_list{e};
    end
end
end

function [correlation, pvalue] = computeCorrelationAndPvalue(counts1, counts2, numPermutations)
correlation = corr( counts1, counts2 );
corrperm = zeros(numPermutations,1);
for i=1:numPermutations; corrperm(i) = corr(counts1, counts2(randperm(numel(counts2)))); end;
%onesidedp = mean(corrperm>c);
%twosidedp = mean(abs(corrperm)>c);
pvalue = mean(abs(corrperm)>correlation);
end