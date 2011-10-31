function prob = spotProb_1D( rna )
% Function for computing the probability that a cell has N number of mRNAs (returned in vector prob) given the probability for each spot to be an mRNA in vector (rna)
%  prob = spotProb( 0.5+0.5*rand(1,10) )

if nargin < 1, rna =  [ .9, .8, .7, .88, .65, .75 ]; end 


Start_Index = 1;
%Initialize the vector of probabilities for different number of mRNAs: 
%{0, 1, 2, .. N}
prob = zeros( 1, length(rna)+1 );

% If the number of spots is large, computing the full Multinomial becomes rather expensive. 
% To avoid combinatoril increase in complexity, only the dimmest N spots
% are considered to be noise (the probabilities for that are computted)
N = 8;
if numel( rna ) > N
	%rna1 = rna;
	% One approach is to use k-mean clusters which is very natural but since the size of each cluster 
	% can be different that requires tracking indecies which is rather complicated 
	%[IDX, rna] = kmeans(rna, N);
	
	%Another, (simpler) approximation is to make the clusters the same size
	%ClusterSize = floor( numel( rna ) / N );
	%Mod =  numel( rna ) - N*ClusterSize;
	%Start_Index = 2+ Mod;
	%if ClusterSize>1
	%	rna = zeros( N+Mod, 1 );	
	%end
	
	%The simplest approach is to assume that the spots with the highest intensities are mRNAs and apply the algorithm 
	% only to the N dimmest spots above a threshold
	vals = sort( rna, 'ascend' );
	Start_Index = numel( rna )-N+1; 
	rna = vals( 1:N );
end

% Vector of probabilities that the spots are not mRNAs
no_rna = 1 - rna;


allInds = 1:length(rna); 

%Probability that the cell has zero RNAs
prob(1) = prod( no_rna );


% Adds path to the combinator function used in the loop below 
%addpath combinator
%addpath bin/combinator

%Probabilities that the cell has i RNAs
for i = 1:length(rna) 	%max(2, length(rna)-10)
	
	 ind = combinator( length(rna), i, 'c' );
	 
	 for j = 1:size( ind, 1)

			prob(Start_Index+i) = ...
			prob(Start_Index+i) +  prod( rna( ind(j,:) ) ) * prod( no_rna(  setdiff( allInds, ind(j,:) ) ) );
	end		


end


% Visualizes the result and checks for accuracy. It should be commented when the function is used for work beyond the development stage 
%plot( prob )
%fprintf( '%1.5f\n',  sum( prob ) )









