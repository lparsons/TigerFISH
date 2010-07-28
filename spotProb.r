#function( rna ) {

rna = c ( .9, .8, .7, .88, .65, .75 );

no_rna = 1 - rna;

#Initialize the vector of probabilities for different number of mRNAs
prob = rep( 0, length(rna)+1 );
allInds = 1:length(rna); 

#Zero RNAs
prob[1] = prod( no_rna );

library( 'gregmisc' )
# 1 RNA
for (i in 2:length(rna)){
	
	 ind = combinations( length(rna), i ); 
	 
	 for (j in 1:attr( ind, "dim")[1] ){

			prob[1+i] = prob[1+i] +  prod( rna[ ind[j,] ] ) * prod( no_rna[  setdiff( allInds, ind[j,] ) ] );
	}		


}

plot( prob, type="b", col = "red",  )
print(  sum( prob ) )

#}
#<environment: namespace:base>	