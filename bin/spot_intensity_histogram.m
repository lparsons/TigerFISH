function f = spot_intensity_histogram(varargin)
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------

ip = inputParser;
ip.FunctionName = 'determine_thresholds';
ip.addRequired('out_spot_intensities',@isnumeric);
ip.addRequired('in_spot_intensities',@isnumeric);
ip.addRequired('threshold',@isnumeric);
ip.addParamValue('max',NaN,@isnumeric);
ip.parse(varargin{:});

% Determine max value and split into 30 bins
%mx = max(vertcat(ip.Results.out_spot_intensities, ip.Results.in_spot_intensities), ip.Results.threshold*6);
% if isnan(ip.Results.max)
    % mx = max( 4*median(ip.Results.in_spot_intensities), ip.Results.threshold*3);
% else
    % mx = ip.Results.max;
% end

mx = max( 3*median(ip.Results.in_spot_intensities(~isnan(ip.Results.in_spot_intensities))),...
	      2*ip.Results.threshold 	);
% if numel(mx) ~= 1 || mx(1) < 0.1
 % mx = 0.1;
% end 

interval = mx/30;
bins = [0:interval:mx,Inf];
    
[i_n i_xout] = histc(ip.Results.in_spot_intensities, bins);
[o_n o_xout] = histc(ip.Results.out_spot_intensities, bins);

% Plot as fraction of all spots
num_spots = (sum(i_n)+sum(o_n));
if ~isempty(num_spots)
    i_f = i_n/(sum(i_n)+sum(o_n));
    o_f = o_n/(sum(i_n)+sum(o_n));
else
    i_f = zeros(size(i_n));
    o_f = zeros(size(o_n));
end
f = figure('Visible', 'off');
if ~isempty(i_f)
    plot(bins, i_f, '-or', 'MarkerFaceColor', 'r', 'LineWidth',3, 'MarkerSize',12);
    hold all
end
if ~isempty(o_f)
    plot(bins, o_f, '--xb', 'MarkerFaceColor', 'b', 'LineWidth',3, 'MarkerSize',12);
end
if ~isnan(ip.Results.threshold)
    ylims = get(gca, 'YLim');
    threshold_label = text(ip.Results.threshold, ylims(1), ...
        num2str(ip.Results.threshold), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');
    threshold_label_extent = get(threshold_label,'Extent');
    threshold_start = threshold_label_extent(2)  + threshold_label_extent(4);
    line([ip.Results.threshold, ip.Results.threshold], [threshold_start, ylims(2)], ...
        'Color', [.3 .3 .3], 'LineStyle', '--', 'LineWidth', 2);
end

hold off
xlim( [0 mx*1.05] ); 
xlabel('Spot Intensity','FontSize', 20,'FontWeight', 'Bold')
ylabel('Fraction','FontSize', 20,'FontWeight', 'Bold')
title('Spot Intensity Histogram','FontSize', 22);
set( gca, 'FontSize', 14, 'FontWeight', 'Bold' ); 
lh = legend( ['Inside Cells: ' num2str(sum(i_n))], ['Outside Cells: '  num2str(sum(o_n))] );
set (lh, 'FontSize', 20,'FontWeight', 'Bold');
    
end
