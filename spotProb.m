function prob = spotProb( rna )

if nargin < 1, rna =  [ .9, .8, .7, .88, .65, .75 ]; end 

no_rna = 1 - rna;

%Initialize the vector of probabilities for different number of mRNAs
prob = zeros( 1, length(rna)+1 );
allInds = 1:length(rna); 

%Zero RNAs
prob(1) = prod( no_rna );


% N RNAs
addpath bin/combinator
for i = max(2, length(rna)-10):length(rna)
	
	 ind = combinator( length(rna), i, 'c' );
	 
	 for j = 1:size( ind, 1)

			prob(1+i) = prob(1+i) +  prod( rna( ind(j,:) ) ) * prod( no_rna(  setdiff( allInds, ind(j,:) ) ) );
	end		


end

plot( prob )
fprintf( '%1.5f\n',  sum( prob ) )
