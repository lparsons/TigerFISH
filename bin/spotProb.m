function prob = spotProb( rna_j, rna_k )
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------

% Function for computing that a cell has N_j mRNAs from the j^th gene and N_k mRNAs from the k^th gene; 


prob_j = spotProb_1D( rna_j );
prob_k = spotProb_1D( rna_k );

prob = prob_j(:) * prob_k; 






