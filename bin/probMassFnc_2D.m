function Y = probMassFnc_2D( x, y, Folder, File_Name, r )  %, Iter
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
    bin_num = numel(bins);

    p_1 = histc(x, bins) * (1/n);  
    p_2 = histc(y, bins) * (1/n);
    
    H_1 = nansum(p_1.*log2(p_1));
    H_2 = nansum(p_2.*log2(p_2));
%end
Y.p_1 = p_1; % marginal dostribution for x
Y.p_2 = p_2; % marginal dostribution for y


% Computes the Probability Mass Function (PMF)
p_1_2 = zeros( bin_num );
% If we have a lot of data bin_num << n so using  bin_num shortens the loop
bins(end) = inf;
for i=1:bin_num-1;        
    p_1_2(:,i) = histc( x( y==bins(i) ),  bins);
end
    p_1_2(:,i+1) = histc( x( y>bins(i) ),  bins);
Y.p_1_2 = p_1_2 * (1/n);


% Plotting the 
sprt = round( MAX/10 );
close all
imagesc( [0 MAX], [MAX 0], log2(p_1_2+1) );
set(gca, 'Ytick', 0:sprt:MAX, 'YtickLabel', MAX:-sprt:0 );
set(gca, 'Xtick', 0:sprt:MAX, 'XtickLabel', 0:sprt:MAX )
%colormap( bone );
BONE = bone;
colormap( BONE(end:-1:1,:) );
 
hc = colorbar;
set(hc, 'Location', 'EastOutside' );
p_1_2_sorted = sort( log2(p_1_2(p_1_2>0)) );
set( gca, 'Clim', [0 log2(p_1_2_sorted(end-2))] );


Segments = regexp(File_Name, '_', 'split' );
%Title = regexprep(File_Name, '_', ' : ' );
% Title = sprintf( '%s : %s :: $R=%+1.2f$',...
%           Segments{2}, ...
%           Segments{3},  r );
Title = sprintf( 'Correlation: $%+1.2f$', r );
          
set( gca, 'FontWeight', 'Bold', 'FontSize', 14 );
h(1) = title( Title ); 
if numel( Segments ) >= 4
h(2) = xlabel( Segments(2) );    
h(3) = ylabel( Segments(3) );
end
set( h(1), 'FontSize', 20 ); % 'interpreter',   'latex' );                                
set(h(2:3), 'fontsize', 18, 'FontWeight', 'bold');

 
% set( gcf, 'Position', [440  359  692   419] )
% set( gca, 'Position', [0.14 0.14 0.80  0.78] );
set( gcf, 'PaperSize', [6.5 6.5],...
          'PaperPositionMode', 'auto' );
print( '-dpdf', [Folder File_Name '.pdf'] );


% Computting Entropy & Mutual Information
H_1_2 = nansum(nansum(Y.p_1_2.*log2(Y.p_1_2)));
Y.MI = H_1_2 - ( H_1 + H_2 );       














