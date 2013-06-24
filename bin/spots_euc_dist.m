function [val ind sz] = spots_euc_dist( x, y )
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------


sz = size(x, 1);
val = zeros(sz,1);
ind = zeros(sz,1);

for i=1:sz
    
   euc_dist = sqrt( ( x(i,1) - y(:,1) ).^2 +...
                    ( x(i,2) - y(:,2) ).^2 +...
                    ( x(i,3) - y(:,3) ).^2      );  
    
   [val(i) ind(i)] = min( euc_dist );  
end