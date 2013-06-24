function sett( Figure_Handle, Font_Size )
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------

if nargin < 2, Font_Size = 18; end 

set( Figure_Handle,  'fontsize',      Font_Size,...
                     'interpreter',   'latex'       );    