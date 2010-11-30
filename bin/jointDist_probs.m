function jointDist_probs( p_1_2,  Gene1, Gene2, Folder_File_Name)  %, Iter
%
% 2) Plots a 2D matrix plot of the PMF and prints it as pdf file


% Plotting the 
MAXy = min( 40, size(p_1_2,2) );  IN.y = 1:MAXy;
MAXx = min( 40, size(p_1_2,1) );  IN.x = 1:MAXx;
spx = round( MAXx/10 );
spy = round( MAXy/10 );

main_fig = figure('Visible', 'off');

sz = 0.65;
sz2   = 0.22;   


h1 = axes( 'Position',    [0.1        0.1       sz         sz ]  );
imagesc( [0 MAXx], [MAXy 0], log2(p_1_2(IN.x,IN.y)+1) ); 

set(h1, 'Ytick', 0:spy:MAXy, 'YtickLabel', MAXy:-spy:0 );
set(h1, 'Xtick', 0:spx:MAXx, 'XtickLabel', 0:spx:MAXx )

%colormap( bone );
BONE = bone;
colormap( BONE(end:-1:1,:) );
 
h_cb = colorbar;
set( h_cb, 'Position', [ sz+0.15   sz+0.12    0.4*sz2    sz2] )
p_1_2_sorted = sort( log2(p_1_2(p_1_2>0)) );
set(h1, 'Clim', [0  1] ); 


%Segments = regexp(File_Name, '_', 'split' );
%Title = regexprep(File_Name, '_', ' : ' );
% Title = sprintf( '%s : %s :: $R=%+1.2f$',...
%           Segments{2}, ...
%           Segments{3},  r );
%Title = sprintf( 'Correlation: $%+1.2f$', r );
          
set(h1, 'FontWeight', 'Bold', 'FontSize', 12 );

h(1) = xlabel( Gene1 );    
h(2) = ylabel( Gene2 );   
%sett( h(1:2), 20 ); 
%set(h(1:2), 'fontsize', 18, 'FontWeight', 'bold');

%h(3) = title( Title ); 
%set( h(3), 'FontSize', 20, 'interpreter',   'latex' );

% bins = 0:MAX;
% Errors1 = sqrt(p1);
% Errors2 = sqrt(p2);
% if 1
%     p1 = log2(p1);   Errors1 = 0.3*p1/2;
%     p2 = log2(p2);   Errors2 = 0.3*p2/2;
% end

p1 = sum( p_1_2, 1 );
p2 = sum( p_1_2, 2 );

h2 = axes( 'Position',    [ 0.1       sz+0.11   sz     sz2 ]  );
bar(1:MAXx, p2(IN.x) ); hold on
%errorbar( 0:MAXx, p2(IN.x), Errors2(IN.x), 'r.' );
ylim( [0 1.02*max( p2(IN.x) )]);
xlim( [1-0.5  MAXx+0.5]);
set(h2, 'Xtick', [] );
set(h2, 'FontWeight', 'Bold' );
h(3) = ylabel( 'Density' );
%==========================================================================




h3 = axes( 'Position',    [sz+0.11    0.1       sz2     sz ]  );
barh(1:MAXy, p1(IN.y) ); hold on
%errorbar_x( p1(IN.y), 0:MAXy, Errors1(IN.y), 'r.' );
xlim( [0 1.02*max( p1(IN.y) )]);
ylim( [1-0.5  MAXy+0.5]);
set(h3, 'Ytick', [] );
set(h3, 'FontWeight', 'Bold' );

h(4) = xlabel( 'Density' ); 
set(h, 'fontsize',   20, 'fontWeight', 'Bold' );
%==========================================================================
%sett( h, 23 ); 


set( gcf, 'Position', [357   245   670   571] ) 
% set( gcf, 'Position', [440  359  692   419] )
% set( gca, 'Position', [0.14 0.14 0.80  0.78] );
set( gcf, 'PaperSize', [7.5  6.5],...
          'PaperPositionMode', 'auto' );
%Hhh = gtext( r  ); sett(Hhh, 34 );    
if nargin >=4
    print(main_fig, '-dpdf', Folder_File_Name );
    %system( [ 'start ' Folder_File_Name '.pdf'] );
end        














