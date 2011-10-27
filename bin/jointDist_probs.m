function Y = jointDist_probs( p_1_2,  Gene1, Gene2, Folder_File_Name)  %, Iter
%
% Plots a 2D matrix plot of the PMF and prints it as pdf file
%
p_1_2 = p_1_2';
% Plotting  
MAXy = min( 40, size(p_1_2,2) );  IN.y = 1:MAXy; Y.y = MAXy;
MAXx = min( 40, size(p_1_2,1) );  IN.x = 1:MAXx; Y.x = MAXx;
spx = round( MAXy/6 );
spy = round( MAXy/6 );

main_fig = figure('Visible', 'on');
set(main_fig,'PaperPositionMode','auto', 'PaperSize', [8  8], 'Units', 'inches')

% Plot main image
h1 = axes( 'Position',    [0.1        0.1       0.65         0.65]  );
imagesc( [0 MAXy], [MAXx 0], log2(p_1_2(IN.x,IN.y)+1) ); 
set(h1, 'Ytick', 0:spx:MAXx, 'YtickLabel', MAXx:-spx:0 );
set(h1, 'Xtick', 0:spy:MAXy, 'XtickLabel', 0:spy:MAXy )
set(h1, 'FontWeight', 'Bold', 'FontSize', 12 );
% Set Clim for imagesc plot
p_1_2_sorted = sort( log2(p_1_2(p_1_2>0)) );
color_lim_max = numel(p_1_2_sorted);
if color_lim_max > 2
    color_lim_max = color_lim_max - 2; % Saturate the top end of the spectrum
end
color_lims = [0 max(1, p_1_2_sorted(color_lim_max))];
set(h1, 'Clim', color_lims ); 
% Set axis lables
h(1) = xlabel( Gene1 );    
h(2) = ylabel( Gene2 );   


% Colormap 
BONE = bone;
colormap( BONE(end:-1:1,:) );
h_cb = colorbar;
set( h_cb, 'Position', [ 0.76   0.76    0.1    0.2] )

         

% Density summary calculation
p1 = sum( p_1_2, 1 );
p2 = sum( p_1_2, 2 );

% Density Summary Y-Axis (Gene2)
h2 = axes( 'Position',    [ 0.1       0.76   0.65     0.2 ]  );
bar(1:MAXy, p1(IN.y) ); hold on
ylim( [0 1.02*max( p1(IN.y) )]);
xlim( [1-0.5  MAXy+0.5]);
set(h2, 'Xtick', [] );
set(h2, 'FontWeight', 'Bold' );
h(3) = ylabel( 'Density' );

% Density Summary X-Axis (Gene1)
h3 = axes( 'Position',    [0.76    0.1       0.2     0.65 ]  );
barh(1:MAXx, p2(IN.x) ); hold on
xlim( [0 1.02*max( p2(IN.x) )]);
ylim( [1-0.5  MAXx+0.5]);
set(h3, 'Ytick', [] );
set(h3, 'FontWeight', 'Bold' );
h(4) = xlabel( 'Density' ); 

% Set fontsize of axis labels
set(h, 'fontsize',   20, 'fontWeight', 'Bold' );

% Print figure if filename specified
if nargin >=4
	set(main_fig,'InvertHardcopy','on')
	set(main_fig, 'PaperPositionMode', 'manual');
	set(main_fig, 'PaperUnits', 'inches');
	set(main_fig, 'PaperPosition', [.25 .25 7.5 7.5]);
    	print(main_fig, '-dpdf', Folder_File_Name, '-r0' );
end        














