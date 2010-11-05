function [cdc_phases cdc_phases_p] = DNA_2_cdc_phases( DNA_Content, MaxIter )



%Sets a Default
if nargin ==1, MaxIter = 1e2; end 


% Generate initial expectations and IC for starting an EM optimization
Median = median( DNA_Content );
muG1 = 0.6*Median;

sigG1 = std( DNA_Content( DNA_Content< 1.4*muG1 ) );



for i=1:MaxIter
    
    %Computes probability for each cells to be in G1 or G2
    
    pG1 = normpdf( DNA_Content,   muG1,   sigG1 ); 
    pG2 = normpdf( DNA_Content, 2*muG1, 2*sigG1 ); 
    
     
    % Uses the probabilities to reassign cells into G1 and G2
    indG1 = find( pG1> pG2 & pG1 > 0.8 );
    indG2 = find( pG1< pG2 & pG2 > 0.8 );
    
    % Updates the expectations muG1 and sigG1 based on 
    muG1 = mean( DNA_Content(indG1) ); 
    muG2 = mean( DNA_Content(indG2) ); 
    muG1 = mean( muG1 + 0.5*muG2 );
    if muG1>Median, muG1=muG1/2; end 
    
    sigG1 = std( DNA_Content(indG1) ); 
    sigG2 = std( DNA_Content(indG2) ); 
    sigG1 = mean( sigG1 + 0.5*sigG2 );
    
    if i>1 && norm(pG1_old-pG1)/norm(pG1) < 1e-2
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

if nargout > 1
    cdc_phases_p.G1 = pG1;
    cdc_phases_p.G2 = pG2;
    cdc_phases_p.muG1 = muG1;
    cdc_phases_p.sigG1 = sigG1;
end 

















