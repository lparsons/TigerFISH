function Y = probMassFnc_2D3( x, y, Folder, File_Name, r )  %, Iter
%
% 1) Computes the marginal distributions of x and y
% as well as thier Probability Mass Function (PMF)
%
% 2) Plots a 2D matrix plot of the PMF and prints it as pdf file
%
% 3) Computes the Mutual Information (MI) between x and y
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------


% If a single argument is passed it has to be nx2 matrix
if nargin < 2, y=x(:,2); x=x(:,1);  end 
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


% Plotting the 
MAXy = min( 40, max(x) );  IN.y = 1:MAXy+1;
MAXx = min( 40, max(y) );  IN.x = 1:MAXx+1;
spx = round( MAXx/10 );
spy = round( MAXy/10 );
close all
sz = 0.65;
sz2   = 0.22;   


h1 = axes( 'Position',    [0.1        0.1       sz         sz ]  );
imagesc( [0 MAXx], [MAXy 0], log2(p_1_2(IN.y,IN.x)+1) ); 

set(h1, 'Ytick', 0:spy:MAXy, 'YtickLabel', MAXy:-spy:0 );
set(h1, 'Xtick', 0:spx:MAXx, 'XtickLabel', 0:spx:MAXx )

%colormap( bone );
BONE = bone;
colormap( BONE(end:-1:1,:) );
 
h_cb = colorbar;
set( h_cb, 'Position', [ sz+0.15   sz+0.12    0.4*sz2    sz2] )
p_1_2_sorted = sort( log2(p_1_2(p_1_2>0)) );
set(h1, 'Clim', [0  p_1_2_sorted(end-4)] ); 


Segments = regexp(File_Name, '_', 'split' );
%Title = regexprep(File_Name, '_', ' : ' );
% Title = sprintf( '%s : %s :: $R=%+1.2f$',...
%           Segments{2}, ...
%           Segments{3},  r );
Title = sprintf( 'Correlation: $%+1.2f$', r );
          
set(h1, 'FontWeight', 'Bold', 'FontSize', 12 );
if numel( Segments ) >= 4
h(1) = xlabel( Segments(2) );    
h(2) = ylabel( Segments(3) );
end   
%sett( h(1:2), 20 ); 
%set(h(1:2), 'fontsize', 18, 'FontWeight', 'bold');

%h(3) = title( Title ); 
%set( h(3), 'FontSize', 20, 'interpreter',   'latex' );

bins = 0:MAX;
Errors1 = sqrt(p1);
Errors2 = sqrt(p2);
if 1
    p1 = log2(p1);   Errors1 = 0.3*p1/2;
    p2 = log2(p2);   Errors2 = 0.3*p2/2;
end

h2 = axes( 'Position',    [ 0.1       sz+0.11   sz     sz2 ]  );
bar(0:MAXx, p2(IN.x) ); hold on
errorbar( 0:MAXx, p2(IN.x), Errors2(IN.x), 'r.' );
ylim( [0 1.02*max( p2(IN.x)+Errors2(IN.x) )]);
xlim( [0-0.5  MAXx+0.5]);
set(h2, 'Xtick', [] );
set(h2, 'FontWeight', 'Bold' );
h(3) = ylabel( '$\log_2(Cells)$' );
%==========================================================================




h3 = axes( 'Position',    [sz+0.11    0.1       sz2     sz ]  );
barh(0:MAXy, p1(IN.y) ); hold on
errorbar_x( p1(IN.y), 0:MAXy, Errors1(IN.y), 'r.' );
xlim( [0 1.02*max( p1(IN.y)+Errors1(IN.y) )]);
ylim( [0-0.5  MAXy+0.5]);
set(h3, 'Ytick', [] );
set(h3, 'FontWeight', 'Bold' );

h(4) = xlabel( '$\log_2(Cells)$' ); sett(h)
%==========================================================================
sett( h, 20 ); 


set( gcf, 'Position', [357   245   670   571] ) 
% set( gcf, 'Position', [440  359  692   419] )
% set( gca, 'Position', [0.14 0.14 0.80  0.78] );
set( gcf, 'PaperSize', [7  6.5],...
          'PaperPositionMode', 'auto' );
print( '-dpdf', [Folder File_Name '.pdf'] );


% Computting Entropy & Mutual Information
H_1_2 = nansum(nansum(Y.p_1_2.*log2(Y.p_1_2)));
Y.MI = H_1_2 - ( H_1 + H_2 );       














