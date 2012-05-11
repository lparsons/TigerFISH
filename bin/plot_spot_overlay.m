function spot_image = plot_spot_overlay(varargin)

ip = inputParser;
ip.FunctionName = 'plot_spot_overlay';
ip.addRequired('spot_data',@isnumeric);
ip.addRequired('threshold',@isnumeric);
ip.addRequired('plot_size',@isnumeric);
ip.addOptional('fgcolor',[1 1 1], @isnumeric);
ip.addOptional('bgcolor',[.4 .4 .4], @isnumeric);
ip.addOptional('N',1, @isnumeric);
ip.parse(varargin{:});

N = ip.Results.N;  % Divider for size until the picture has a "good" size on screen

% Threshold spots
xx0 = ip.Results.spot_data(ip.Results.spot_data(:,4)<=ip.Results.threshold,1);
yy0 = ip.Results.spot_data(ip.Results.spot_data(:,4)<=ip.Results.threshold,2);
xx1 = ip.Results.spot_data(ip.Results.spot_data(:,4)>ip.Results.threshold,1);
yy1 = ip.Results.spot_data(ip.Results.spot_data(:,4)>ip.Results.threshold,2);

% xx1 = spot_data(spot_data(:,4)>thresholds(1)...
%     &spot_data(:,4)<=thresholds(2),1);
% yy1 = spot_data(spot_data(:,4)>thresholds(1)...
%     &spot_data(:,4)<=thresholds(2),2);
% xx2 = spot_data(spot_data(:,4)>thresholds(2),1);
% yy2 = spot_data(spot_data(:,4)>thresholds(2),2);

% Generate figure
spot_image = figure('Visible', 'off');
set(spot_image, 'Units', 'pixels', 'Position', [0, 0, ip.Results.plot_size(2) / N, ip.Results.plot_size(1) / N]);
set(gca, 'Units', 'normalized', 'Position', [0,0,1,1]);
%set(gca, 'Units', 'points');

hold on;
plot(xx0,yy0,'o','MarkerSize',10/N, 'Color',ip.Results.bgcolor);
plot(xx1,yy1,'o','MarkerSize',10/N, 'Color',ip.Results.fgcolor);
%plot(xx2,yy2,'o','MarkerSize',4, 'MarkerColor',ip.Results.bgcolor);
set(gca,'YDir','reverse');
xlim([0 ip.Results.plot_size(2)]); ylim([0 ip.Results.plot_size(1)]);
%axis image; 
hold off;

% Set figure options
%set(spot_image, 'Units', 'points', 'PaperUnits', 'points', 'PaperPositionMode', 'auto');
set(spot_image,'Color','black');
set(gca,'Color','black');
set(spot_image,'InvertHardCopy','off');

% Print figure
% screen_DPI = get(0, 'ScreenPixelsPerInch');
% print(spot_image, '-dpng', sprintf('-r%d', N * screen_DPI), [this_exp_reg_results_dir filesep 'region_' num2str(region_number) '_' dye '_spot_image.png']);
% % Load and save image with transparency
% tmp = imread([this_exp_reg_results_dir filesep 'region_' num2str(region_number) '_' dye '_spot_image.png']);
% imwrite(tmp,[this_exp_reg_results_dir filesep 'region_' num2str(region_number) '_' dye '_spot_image.png'], 'png','Transparency', [0,0,0]);
% 
% close(spot_image);

end
