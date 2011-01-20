function Y = jointDist( x, y,  Gene1, Gene2, Folder_File_Name, Ylim)  %, Iter
%
% 1) Computes the marginal distributions of x and y
% as well as thier Probability Mass Function (PMF)
%
% 2) Plots a 2D matrix plot of the PMF and prints it as pdf file
%
% 3) Computes the Mutual Information (MI) between x and y


% If a single argument is passed it has to be nx2 matrix
if   min(size(x)) == 2, y=x(:,2); x=x(:,1);  end 
%If no iteration is passed the function will compute the
% marginal distr. 
%if nargin < 3, Iter = 1; end 

%Makes sure the input vectors have the same size
l1 = length(x);
l2 = length(y);
if l1 ~= l2, error( 'Different Sizes !' ); else n = l1; end 
%% 

% Computes the marginal dostributions.
% We should avoid recomputing the marginal dostributions thousands of times when bootstrapping
%if Iter == 1     
    %bins = unique([x y]);
    MAX = max( max(x), max(y) );
    MAX = min( 40, MAX );
    bins = 0:MAX;
    %bins(end) = inf;
    bin_num = numel(bins);

    p1 = histc(x, bins); 
    p2 = histc(y, bins);
    
    p_1 = p1  * (1/n); 
    p_2 = p2  * (1/n); 
     
    H_1 = nansum(p_1.*log2(p_1));
    H_2 = nansum(p_2.*log2(p_2));
%end
Y.p_1 = p_1; % marginal dostribution for x
Y.p_2 = p_2; % marginal dostribution for y


% Computes the Probability Mass Function (PMF)
p_1_2 = zeros( bin_num );
% If we have a lot of data bin_num << n so using  bin_num shortens the loop
for i=1:bin_num        
    p_1_2(:,i) = histc( x( y==bins(i) ),  bins);
end
    %p_1_2(:,i+1) = histc( x( y>bins(i) ),  bins);
    p1 = sum(p_1_2,2);
    p2 = sum(p_1_2,1);
Y.p_1_2 = p_1_2 * (1/n);


% Plotting
if nargin < 6
    MAXx = min( 40, max(x) );  IN.x = 1:MAXx+1;    
    MAXy = min( 40, max(y) );  IN.y = 1:MAXy+1;
else
    MAXy = Ylim.y;  IN.y = 1:MAXy+1;
    MAXx = Ylim.x;  IN.x = 1:MAXx+1;    
end
spx = round( MAXx/10 );
spy = round( MAXy/10 );

main_fig = figure('Visible', 'on');
%set(main_fig,'PaperPositionMode','auto', 'PaperSize', [8  8], 'Units', 'inches')


% Plot main image
h1 = axes( 'Position',    [0.1       0.1      0.65     0.65 ], 'Parent', main_fig, 'Units', 'normalized'  );
imagesc( [0 MAXx], [MAXy 0], log2(p_1_2(IN.x,IN.y)+1)', 'Parent',h1 ); 
set(h1, 'Ytick', 0:spy:MAXy, 'YtickLabel', MAXy:-spy:0 );
set(h1, 'Xtick', 0:spx:MAXx, 'XtickLabel', 0:spx:MAXx )
set(h1, 'FontWeight', 'Bold', 'FontSize', 12 );
% Set Clim for imagesc plot
%p_1_2
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
set( h_cb, 'Position', [ 0.76   0.76    0.1    0.2] , 'Parent', main_fig, 'Units', 'normalized')
caxis(color_lims)

% Log_2(Cells) Summary calculation
bins = 0:MAX;
% Errors1 = sqrt(p1);
% Errors2 = sqrt(p2);
% if 1
%     p1 = log2(p1);   Errors1 = 0.3*p1/2;
%     p2 = log2(p2);   Errors2 = 0.3*p2/2;
% end

% Log_2(Cells) Summary Y-Axis (Gene2)
h2 = axes( 'Position',    [ 0.1      0.76   0.65     0.2 ] , 'Parent', main_fig, 'Units', 'normalized' );
bar(0:MAXx, p1(IN.x) ); hold on
ylim( [0 1.02*max( p1(IN.x) )]);
xlim( [0-0.5  MAXx+0.5]);
set(h2, 'Xtick', [] );
set(h2, 'FontWeight', 'Bold' );
h(3) = ylabel( 'log_2(Cells)' );


% Log_2(Cells) Summary X-Axis (Gene1)
h3 = axes( 'Position',    [0.76    0.1       0.2     0.65 ] , 'Parent', main_fig, 'Units', 'normalized' );
barh(0:MAXy, p2(IN.y) ); hold on
xlim( [0 1.02*max( p1(IN.y) )]);
ylim( [0-0.5  MAXy+0.5]);
set(h3, 'Ytick', [] );
set(h3, 'FontWeight', 'Bold' );
h(4) = xlabel( 'log_2(Cells)' ); 


% Set fontsize of axis labels
set(h, 'fontsize',   20, 'fontWeight', 'Bold' );


% Print figure if filename specified
if nargin >=5
	set(main_fig,'InvertHardcopy','on')
	set(gcf, 'PaperPositionMode', 'manual');
	set(gcf, 'PaperUnits', 'inches');
	set(gcf, 'PaperPosition', [.25 .25 7.5 7.5]);
	print(main_fig, '-dpdf', Folder_File_Name, '-r0' );
end 

% Computting Entropy & Mutual Information
H_1_2 = nansum(nansum(Y.p_1_2.*log2(Y.p_1_2)));
Y.MI = H_1_2 - ( H_1 + H_2 );       



