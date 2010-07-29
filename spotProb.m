function prob = spotProb( rna )

% Function for computing the probability that a cell has N number of mRNAs (returned in vector prob) given the probability for each spot to be an mRNA in vector (rna)
%  prob = spotProb( 0.5+0.5*rand(1,10) )

if nargin < 1, rna =  [ .9, .8, .7, .88, .65, .75 ]; end 

% If the number of spots is large, computing the full Multinomial becomes rather expensive. 
% To avoid combinatoril increase in complexity, the spots are grouped into N number of grous
N = 12;
if numel( rna ) > N
	rna1 = rna;
	[IDX, rna] = kmeans(rna, N)
end

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
