function [spot_data sp potSpot] = measure_spots_np( spot_locations, img, Thresh  )


if nargin <= 2, Thresh = 1.10; end
[sz.x sz.y sz.z] = size( img.layers );
for i=1:max( size( img.info ) )
    layer = img.layers(:,:,i);
    img.layers(:,:,i) = layer * ( numel(layer) / sum(layer(:)) ); %median(layer(:)); %
end


ind = find( spot_locations );
[y x z] = ind2sub( size(spot_locations), ind );



xyz = [x y z]; 
spot_data = zeros( size(xyz,1), 4 ); 
spot_data(:,1) = x;
spot_data(:,2) = y;
spot_data(:,3) = z;


%Removes Border Pixels 
xyz = xyz( x > 6      & y > 6      & z > 3 &...
           x < sz.x-5 & y < sz.y-5 & z < sz.z-2,   :  );

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
i_num5 = 1/( r3*r5^2 - r2*r4^2 ); j=0;
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
    

    
    [sl hat_xyz Rsq] = spot_regress_fun( pix_53, 11 );
    
    

    Mean = sum( pix_3(:) ) * i_num3;
    Background = (sum(pix_6(:)) - sum(pix_5(:)) ) * i_num6;
        
    Contrast = Mean/Background;
    potSpot.contrast(i) = Contrast;
    potSpot.rsq( i ) = Rsq;
    potSpot.sl( i,: ) = sl;
    

%     if 1 && (Rsq > 0.10 && sl(1)<0 ) || ...
%        Contrast > Thresh 
%   
       j=j+1;       
 
       sp.mean(j) =  Mean;      
       sp.background(j) = Background;
       pix = pix_5(:) - Background; %Background_median;
       sp.intensity(j, 1) = ( pix' * pix ) / sum( pix );
       sp.xyz(j,:) = xyz(i,:);
       %sp.pixs{j} = pix_b;
       sp.rsq(j) = Rsq;
       sp.hat_xyz(j,:) = hat_xyz;
       sp.contrast(j) = Contrast;
       
       spot_data(i,4) = sp.intensity(j, 1);
%    end
end
%fprintf( '3) After Considering the Neigbourhoods: \t%d \n', j);

%if j <= 2, sp.contrast = 0; spb.contrast = 0; return; end 


return 



%% Some spots may come from the same location. This module 
% secelts the best spot from each uniqe location  

To_Be_Skipped = [];
bestUniqe = [];
sz = size(sp.xyz,1);
for i=1:sz

   if sum( To_Be_Skipped == i ) == 1
       continue 
   end 

   euc_dist = sqrt( ( sp.xyz(i,1) - sp.xyz(:,1) ).^2 +...
                    ( sp.xyz(i,2) - sp.xyz(:,2) ).^2 +...
                    ( sp.xyz(i,3) - sp.xyz(:,3) ).^2    );  

   ind = find( euc_dist < 8 );

   [val in] = max( sp.contrast(ind) ); 

   in = ind(in); 

   To_Be_Skipped =[To_Be_Skipped;  ind(ind~=in) ];

   if sum( bestUniqe == in ) == 0
      bestUniqe = [bestUniqe; in];
   end
end
sp.bestUniqe = bestUniqe;
fprintf( '4) After Filtering Non Uniqe Spots: \t%d \n', numel(bestUniqe));
if numel( bestUniqe ) > 1500
    Intensities = sp.intensity( bestUniqe ); 
    [sorted_intensity, ...
     sorted_intensity_vals ] = sort( Intensities, 'descend' ); 

    bestUniqe = bestUniqe( sorted_intensity_vals(1:1000) );
end 
spb.intensity = sp.intensity( bestUniqe );
spb.contrast = sp.contrast( bestUniqe );
spb.rsq = sp.rsq( bestUniqe );
spb.xyz = sp.xyz( bestUniqe, : );


return

if numel(bestUniqe) <= 2, sp.contrast = 0; spb.contrast = 0; return; end 
%%
[spb.km.inds,...
 spb.km.means,...
 spb.km.sums  ] =  kmeans(spb.intensity(:), 2, 'replicates', 4);

[spb.km2.inds,...
 spb.km2.means,...
 spb.km2.sums  ] =  kmeans( [spb.intensity(:)...
                             spb.contrast(:) ], 2, 'replicates', 4);
                         
[val ind] = max( spb.km.means );
if  diff( spb.km.means ) / min(spb.km.means) >2 || ...
    sum( spb.km.inds == ind ) > sum( spb.contrast > 1.20 )

         spb.km.mRNAs_ind = find( spb.km.inds == ind ); 
else
         %spb.km.mRNAs_ind = find( spb.intensity(:) >= 1.20  ); 
         spb.km.mRNAs_ind = 1:numel(bestUniqe);
end 

spb.mat = [ spb.xyz( spb.km.mRNAs_ind, : )  spb.intensity( spb.km.mRNAs_ind ) ];

spb.median = median( img.layers(:) );
return