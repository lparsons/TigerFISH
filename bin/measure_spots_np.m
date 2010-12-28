function [spot_data sp potSpot] = measure_spots_np( spot_locations, img, Thresh  )


if nargin <= 2, Thresh = 1.10; end
[sz.x sz.y sz.z] = size( img.layers );
for i=1:max( size( img.info ) )
    layer = img.layers(:,:,i);
    img.layers(:,:,i) = layer * ( numel(layer) / sum(layer(:)) ); %median(layer(:)); %
end


ind = find( spot_locations );
[x y z] = ind2sub( size(spot_locations), ind );



%Removes Border Pixels 
xyz = [x y z];
xyz = xyz( x > 6      & y > 6      & z > 3 &...
           x < sz.x-5 & y < sz.y-5 & z < sz.z-2,   :  );

spot_data = zeros( size(xyz,1), 4 ); 
spot_data(:,1) = xyz(:,2);
spot_data(:,2) = xyz(:,1);
spot_data(:,3) = xyz(:,3);
       
Num =  size( xyz, 1 );    

if Num==0, sp = []; return, end

potSpot.contrast = zeros( Num, 1 );
potSpot.size = zeros( Num, 1 );
potSpot.rsq = zeros( Num, 1 );
potSpot.sl = zeros( Num, 2 );

rg1 = -1:1;   r1 = 2*1+1;
rg2 = -2:2;   r2 = 2*2+1;
rg3 = -3:3;   r3 = 2*3+1;
rg4 = -4:4;   r4 = 2*4+1;
rg5 = -5:5;   r5 = 2*5+1;
rg6 = -6:6;   r6 = 2*6+1;

i_num3 = 1/( r2*r3^2 ); % volume in pixels 
i_num4 = 1/( r2*r4^2 ); % volume in pixels 
i_num5 = 1/( r3*r5^2 - r2*r4^2 );
i_num6 = 1/( r3*r6^2 - r2*r5^2 );
%Background_median = median( img.layers(:) );
for i=1:Num
    
    x = xyz(i,1);
    y = xyz(i,2);
    z = xyz(i,3);

	%pix_1 = img.layers( x+rg1, y+rg1, z+rg1 ); % cube
	%pix_2 = img.layers( x+rg2, y+rg2, z+rg2 ); % cube
	pix_3 = img.layers( x+rg3, y+rg3, z+rg2 ); % cube
    pix_4 = img.layers( x+rg4, y+rg4, z+rg2 );
	pix_5 = img.layers( x+rg5, y+rg5, z+rg2 );
    pix_53 = img.layers( x+rg5, y+rg5, z+rg3 );
    pix_6 = img.layers( x+rg6, y+rg6, z+rg3 );
    
    
%     if xyz(i,3) > 4 && xyz(i,3) < sz.z-4 
% 
%       pix_5 = img.layers( x+rg5, y+rg5, z+rg3 );
% 
%       Z = z+rg4;  i_num_ = i_num5;  Sz = 11;
%     else
%       X = x+rg4;
%       Y = y+rg4;
%       Z = z+rg3;  i_num_ = i_num4;  Sz = 9;
%     end
    

    
    %[sl hat_xyz Rsq] = spot_regress_fun( pix_53, 11 );
    
    

    %Mean = sum( pix_3(:) ) * i_num3;
    Background = (sum(pix_6(:)) - sum(pix_5(:)) ) * i_num6;
        
    %Contrast = Mean/Background;
    %potSpot.contrast(i) = Contrast;
    %potSpot.rsq( i ) = Rsq;
    %potSpot.sl( i,: ) = sl;
    
 
    %sp.mean(i) =  Mean;
    sp.background(i) = Background;
    pix = pix_5(:) - Background; %Background_median;
    sp.intensity(i, 1) = ( pix' * pix ) / sum( pix );
    %sp.xyz(i,:) = xyz(i,:);
    %sp.pixs{j} = pix_b;
    %sp.rsq(i) = Rsq;
    %sp.hat_xyz(i,:) = hat_xyz;
    %sp.contrast(i) = Contrast;
    
    spot_data(i,4) = sp.intensity(i, 1);
end


return 