function spot_locations = filter_border_spots( spot_locations )


[sz.x sz.y sz.z] = size( spot_locations );
ind = find( spot_locations );
[x y z] = ind2sub( size(spot_locations), ind );

xyz = [x y z]; 

%Removes Border Pixels 
Indx = find( x > 6      & y > 6      & z > 3 &...
             x < sz.x-5 & y < sz.y-5 & z < sz.z-2  );

INDEX = sub2ind( size(spot_locations), x(Indx),  y(Indx),  z(Indx)  );        

spot_locations = zeros( size(spot_locations) );

spot_locations(INDEX) = 1; 