% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------
%CUMSUMALL cumulative sum of integer elements 
%   For vectors, CUMSUMALL(X) is a vector containing the cumulative sum of
%   the elements of X. For matrices, CUMSUM(X) is a matrix the same size
%   as X containing the cumulative sums over each column.  
% 
% Class suport:
%     int8, int16, int32
%
% Keep in mind that the usefullness of this MEX-File is limited because of
% saturation for most problems outside of use with COMBINATOR.
%
%   See also cumsum, CUMPROD, SUM, PROD.
