function prob = spotProb( rna_j, rna_k )

% Function for computing that a cell has N_j mRNAs from the j^th gene and N_k mRNAs from the k^th gene; 


prob_j = spotProb_1D( rna_j );
prob_k = spotProb_1D( rna_k );

prob = prob_j(:) * prob_k; 






