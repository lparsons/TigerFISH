function f = spot_intensity_histogram(varargin)

ip = inputParser;
ip.FunctionName = 'determine_thresholds';
ip.addRequired('out_spot_intensities',@isnumeric);
ip.addRequired('in_spot_intensities',@isnumeric);
ip.addRequired('threshold',@isnumeric);
ip.addParamValue('max',NaN,@isnumeric);
ip.parse(varargin{:});

% Determine max value and split into 30 bins
%mx = max(vertcat(ip.Results.out_spot_intensities, ip.Results.in_spot_intensities), ip.Results.threshold*6);
if isnan(ip.Results.max)
    mx = max( 4*median(ip.Results.in_spot_intensities), ip.Results.threshold*3);
else
    mx = ip.Results.max;
end
interval = mx/30;
bins = [0:interval:mx,Inf];
    
[i_n i_xout] = histc(ip.Results.in_spot_intensities, bins);
i_f = i_n/sum(i_n);
[o_n o_xout] = histc(ip.Results.out_spot_intensities, bins);
o_f = o_n/sum(o_n);

f = figure('Visible', 'off');
plot(bins, i_f, '-or', 'MarkerFaceColor', 'r', 'LineWidth',3, 'MarkerSize',12);
hold all
plot(bins, o_f, '--xb', 'MarkerFaceColor', 'b', 'LineWidth',3, 'MarkerSize',12);
if ~isnan(ip.Results.threshold)
    ylims = get(gca, 'YLim');
    line([ip.Results.threshold, ip.Results.threshold], ylims, 'Color', [.3 .3 .3], 'LineStyle', '--', 'LineWidth', 2);
end

hold off
xlabel('Spot Intensity','FontSize', 20,'FontWeight', 'Bold')
ylabel('Fraction','FontSize', 20,'FontWeight', 'Bold')
title('Spot Intensity Histogram','FontSize', 22);
set( gca, 'FontSize', 14, 'FontWeight', 'Bold' ); 
lh = legend('Inside Cells', 'Outside Cells');
set (lh, 'FontSize', 20,'FontWeight', 'Bold');
    
end