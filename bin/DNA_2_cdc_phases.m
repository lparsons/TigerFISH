function [cdc_phases cdc_phases_p] = DNA_2_cdc_phases( DNA_Content, MaxIter, PathFileName )

% Set MaxIter default
if nargin ==1 || isempty(MaxIter), MaxIter = 1e2; end 

% Normalize variance to avoid problems from treating pdf as pmf
Std = std(  DNA_Content( ~isnan(DNA_Content) )  );
DNA_Content = DNA_Content * (1/Std );
var(DNA_Content);
 
% Generate initial expectations and IC for starting an EM optimization
Median = median( DNA_Content );
muG1 = 0.6*Median;
sigG1 = std( DNA_Content( DNA_Content< 1.4*muG1 ) );

for i=1:MaxIter
    
    % Computes probability for each cells to be in G1 or G2
    pG1 = normpdf( DNA_Content,   muG1,   sigG1 ); 
    pG2 = normpdf( DNA_Content, 2*muG1, 2*sigG1 ); 
    
    % Computes cumulative probability for each cells to be in G1 or G2
    pcG1 = normcdf( DNA_Content,   muG1,   sigG1 ); 
    pcG2 = normcdf( DNA_Content, 2*muG1, 2*sigG1 );     
     
    % Uses the probabilities to reassign cells into G1 and G2
    indG1b = find( DNA_Content< muG1 & pcG1 > 0.01 );    
    indG1a = find( DNA_Content> muG1 & pcG1 < 0.8  & pG1> pG2 );
    indG1 = [indG1a; indG1b];
    indG2b = find( DNA_Content> 2*muG1 & pcG2 < 0.99 );    
    indG2a = find( DNA_Content< 2*muG1 & pcG2 > 0.2  & pG1< pG2 );
    indG2 = [indG2a; indG2b];
    
    % Updates the expectations muG1 and sigG1 based on 
    muG1 = mean( DNA_Content(indG1) ); 
    muG2 = mean( DNA_Content(indG2) ); 
    muG1 = mean( muG1 + 0.5*muG2 );        % fprintf( '%1.2f\n', muG1 );
    if muG1>Median, muG1=muG1/2; end
    if muG1<0.1*Median, muG1=2*muG12; end
    
    sigG1 = std( DNA_Content(indG1) ); 
    sigG2 = std( DNA_Content(indG2) ); 
    sigG1 = mean( sigG1 + 0.5*sigG2 );
    
    if i>1 && norm(pG1_old-pG1)/norm(pG1) < 1e-3
        break
    end
    pG1_old = pG1;
end 

cdc_phases = zeros( size(DNA_Content) );  
cdc_phases( indG1 ) = 1;
cdc_phases( indG2 ) = 3;

indS = find( cdc_phases == 0       &...
             DNA_Content >    muG1 &...
             DNA_Content <  2*muG1      );
         
cdc_phases( indS ) = 2;


%Returns probabilities if requested
if nargout > 1
    cdc_phases_p.G1 = pG1;
    cdc_phases_p.G2 = pG2;
    cdc_phases_p.muG1 = muG1;
    cdc_phases_p.sigG1 = sigG1;
    cdc_phases_p.DNA_Content = DNA_Content;
end 


%Prints distribution of DNA content and separation by phases if filename
%(and path a provided as a third argument)

if nargin >= 3 && ~isempty(PathFileName)

    figure('Visible', 'off')
    
    l = linspace( min(0,min(DNA_Content)),...
                 min(3*Median, max(DNA_Content)), 30 ); clear fr  

    frBar = histc( DNA_Content, l );
    hbar = bar( l, frBar ); hold on
    set(hbar,...
        'FaceColor', 'w',...
        'EdgeColor', 'k',...
        'BarWidth', 0.8  );
    for i=1:3
      fr(:,i) = histc( DNA_Content(cdc_phases==i), l );
    end

    plot( l, fr, 'LineWidth', 3 )
    fprintf('CDC Phases xlim %f : %f', l(1), l(end));
    fprintf(  '%1.2f\t', DNA_Content);
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

    print( '-dpdf', PathFileName );
      
end

