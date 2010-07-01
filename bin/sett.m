function sett( Figure_Handle, Font_Size )

if nargin < 2, Font_Size = 18; end 

set( Figure_Handle,  'fontsize',      Font_Size,...
                     'interpreter',   'latex'       );    