function [mat spb sp img] = locSpots_2( X )
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------

%X = '/data/spots/E191_R3_ssk22-16,95,250,368_ssk22-1_yef3/Image_Cy3_001.tif';

if isnumeric(X),
   if namel(size(X)) ~= 3, 
       error( 'Incorrect Input to locSpots_2 !' ); 
   end
   for i=1:size(X,3) 
     img.layers(:,:,i) = X(:,:,i) * (1/sum(sum(X(:,:,i))));
   end
elseif ischar(X),  filename = X;     % Read img layers
    img.info = imfinfo( filename );
    img.layers = zeros(img.info(1).Height,...
                         img.info(1).Width,...
                         size(img.info, 1)    );
    for i=1:max( size( img.info ) )
        layer = double( imread( filename, i ) );
        img.layers(:,:,i) = layer * ( numel(layer) / sum(layer(:)) );
    end
    
end 
[sz.x sz.y sz.z] = size( img.layers );
%%
[img.Max.Vals  img.Max.Inds] = max( img.layers, [], 3 );
 Var = var( img.layers, [], 3 );

[x y] = find( Var>2.5*mean(Var(:))  ); %& Max.Vals > 1
if numel(x) > 2e4 ||  numel(x) < 3e3
    [sorted_vals sorted_inds]= sort( Var(:), 'descend' );
    [x y] = ind2sub( size(Var), sorted_inds(1:8e3) ); 
end
fprintf( 'Number of Selected Locations: \n' )
fprintf( '1) High Varaince: \t%d \n',  size(x,1) );   %Intensity and

xy = sub2ind( [sz.x sz.y], x, y );
z = img.Max.Inds(xy);
xyz = [x y z];

%Removes Border Pixels 
xyz = xyz( x > 3      & y > 3      & z > 3 &...
           x < sz.x-2 & y < sz.y-2 & z < sz.z-2,   :  );
Num = size(xyz,1);
fprintf( '2) After Removal of Border Pixels: \t%d \n',  Num );
%% 
potSpot.contrast = zeros( Num, 1 );
potSpot.size = zeros( Num, 1 );
potSpot.rsq = zeros( Num, 1 );
potSpot.sl = zeros( Num, 2 );

r1 = 1;  rg1 = -r1:r1;   rr1 = 2*r1+1;
r2 = 2;  rg2 = -r2:r2;   rr2 = 2*r2+1;
r3 = 3;  rg3 = -r3:r3;   rr3 = 2*r3+1;
r4 = 4;  rg4 = -r4:r4;   rr4 = 2*r4+1;

i_num1 = 1/( rr1^3 );
i_num3 = 1/( rr3^3 - rr2^3 ); 
i_num4 = 1/( rr4^3 - rr2^3 ); j=0;
for i=1:Num
    
    x = xyz(i,1);
    y = xyz(i,2);
    z = xyz(i,3);
    
    if sum( xyz(i,:)-4 >= ones(1,3) ) == 3 && ...
       sum( xyz(i,:)+4 <= [sz.x sz.y sz.z] ) == 3  
      X = x+rg4;
      Y = y+rg4;
      Z = z+rg4;  i_num_ = i_num4;  Sz = 9;
    else
      X = x+rg3;
      Y = y+rg3;
      Z = z+rg3;  i_num_ = i_num3;  Sz = 7;
    end
    
    pix_1 = img.layers( x+rg1, y+rg1, z+rg1 );
    pix_2 = img.layers( x+rg2, y+rg2, z+rg2 );
    pix_b = img.layers( X, Y, Z );
    
    [sl hat_xyz Rsq] = spot_regress_fun( pix_b, Sz );
    
    

    Mean = sum( pix_1(:) ) * i_num1;
    Background = (sum(pix_b(:)) - sum(pix_2(:)) ) * i_num_;
        
    Contrast = Mean/Background;
    potSpot.contrast = Contrast;
    potSpot.size( i ) = Sz;
    potSpot.rsq( i ) = Rsq;
    potSpot.sl( i,: ) = sl;
    

    if Rsq      > 0.35 || ...
       Contrast > 1.15 
  
       j=j+1;       
 
       sp.mean(j) =  Mean;      
       sp.background(j) = Background;
       pix = pix_2(:) - Background;
       sp.intensity(j) = ( pix' * pix ) / sum( pix );
       sp.size(j) = Sz;
       sp.xyz(j,:) = xyz(i,:);
      %sp.pixs{j} = pix_b;
       sp.rsq(j) = Rsq;
       sp.hat_xyz(j,:) = hat_xyz;
       sp.contrast(j) = Contrast;
   end
end
fprintf( '3) After Considering the Neigbourhoods: \t%d \n', j);

if j <= 2, sp.contrast = 0; spb.contrast = 0; return; end 

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
if numel( bestUniqe ) > 500
    Intensities = sp.intensity( bestUniqe ); 
    [sorted_intensity, ...
     sorted_intensity_vals ] = sort( Intensities, 'descend' ); 

    bestUniqe = bestUniqe( sorted_intensity_vals(1:400) );
end 
spb.intensity = sp.intensity( bestUniqe );
spb.contrast = sp.contrast( bestUniqe );
spb.rsq = sp.rsq( bestUniqe );
spb.xyz = sp.xyz( bestUniqe, : );

fprintf( '4) After Filtering Non Uniqe Spots: \t%d \n', numel(bestUniqe));

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

mat = [ spb.xyz( spb.km.mRNAs_ind, : )  spb.intensity( spb.km.mRNAs_ind ) ];


return
%% == For Interactive Exploration of the Results ==
%%
ind = find( potSpot.contrast > 1.1 | ...
            potSpot.rsq      > 0.3     );

fprintf( '3) After Considering the Neigbourhoods: \t%d \n', ...
          numel(ind)    );
%%      
ii=0;
ind = find( potSpot.should_Be_Counted );
%ind = find( potSpot.contrast<1.12 & potSpot.rsq>0.38 );
%ind = find( potSpot.rsq < 0.2 & potSpot.contrast>1.2);
%ind = find( potSpot.rsq > 0.3 & potSpot.std>0.2)
%%   Visualization 
ii=ii+1;                                             

i = ind(ii);
x = potSpot.xyz(i,1)+rg4;
y = potSpot.xyz(i,2)+rg4;
z = potSpot.xyz(i,3)+rg4;
potSpot.xyz(i,3), pause( 1 )

pix = img.layers( x, y, z );

%[Mean Cov Rsq] = mle_3d_Gaussian(   pix  )

sorted_pix = sort( pix(:) );
pix = pix(:) - sorted_pix(1);
intensity = ( pix' * pix ) / sum( pix )

Clim = [sorted_pix(50), sorted_pix(end-50)];

fig
for j = z(1):z(end)
    imagesc( img.layers( x, y, j ) ); colorbar

    set( gca, 'Clim', Clim );
    colormap( bone ); %redgreencmap
    colorbar
    set( gcf, 'Name', sprintf( 'Layer: %d', j ) ); 
    pause( 2 )
end  
close all  

%%





















