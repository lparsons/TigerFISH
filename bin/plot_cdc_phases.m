function [ f ] = plot_cdc_phases( dna_content, cdc_phases, filename )
%PLOT_CDC_PHASES Prints distribution of DNA content and phases
%   Prints distribution of DNA content and separation by phases
%   Returns figure handle

% Normalize variance to avoid problems from treating pdf as pmf
dna_content_std = std(  dna_content( ~isnan(dna_content) & ~isinf(dna_content) ) );
dna_content_norm = dna_content * (1/dna_content_std );
dna_content_median = median( dna_content( ~isnan(dna_content) & ~isinf(dna_content) ) );

f = figure('Visible', 'off');

dna_content_nonan = dna_content_norm;
dna_content_nonan( isnan(dna_content_nonan) ) = 0;

l = linspace( min(0,min(dna_content_nonan)),...
    min(3*dna_content_median, max(dna_content_nonan)), 30 ); clear fr

l = unique(l); % Prevent bar plot function from crashing

frBar = histc( dna_content_nonan, l );

% TODO - This often fails and crashed
% Error: XData cannot contain duplicate values)
hbar = bar( l, frBar ); hold on

set(hbar,...
    'FaceColor', 'w',...
    'EdgeColor', 'k',...
    'BarWidth', 0.8  );
for i=1:3
    fr(:,i) = histc( dna_content_nonan(cdc_phases==i), l );
end

plot( l, fr, 'LineWidth', 3 )
%fprintf('CDC Phases xlim %f : %f\n', l(1), l(end));
%fprintf(  '%1.2f\t\n', dna_content);
xlim( [l(1) l(end)] );
h(1) = xlabel( 'DNA Content' );
h(2) = ylabel( 'Number of Cells' );
h(3) = legend( 'All', 'G1', 'S', 'G2'  );
set(h, 'fontsize',       24,...
    'interpreter',   'tex'       );

%set( gcf, 'Position', [440  159  692   419] )
%set( gca, 'Position',  [0.14 0.14 0.80  0.78] );
set( gca, 'FontSize', 14, 'FontWeight', 'Bold' );
set( gcf, 'PaperSize', [8 5],...
    'PaperPositionMode', 'manual',...
    'PaperPosition', [.1 .1 7.8 4.8] );

print( '-dpdf', filename );

end

