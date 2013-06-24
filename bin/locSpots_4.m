function [spb sp img] = locSpots_4( X, Thresh )
% ------------------------
% Copyright (c) 2010-2013, Lance R. Parsons <lparsons@princeton.edu>, Nikolai Slavov <nslavov@mit.edu>
% Licensed under the BSD 2-Clause License: http://opensource.org/licenses/BSD-2-Clause
% ------------------------

%X = '/data/spots/E191_R3_ssk22-16,95,250,368_ssk22-1_yef3/Image_Cy3_001.tif';
if nargin <= 1, Thresh = 1.1; end 

if isnumeric(X),
   if numel(size(X)) ~= 3, 
       error( 'Incorrect Input to locSpots_2 !' ); 
   end
   for i=1:size(X,3) 
     img.layers(:,:,i) = X(:,:,i) * (numel(X(:,:,i))/sum(sum(X(:,:,i))));
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
xyz = xyz( x > 5      & y > 5      & z > 4 &...
           x < sz.x-4 & y < sz.y-4 & z < sz.z-3,   :  );
Num = size(xyz,1);
fprintf( '2) After Removal of Border Pixels: \t%d \n',  Num );
%% 
potSpot.contrast = cell( Num, 1 );
potSpot.size = zeros( Num, 1 );
potSpot.rsq = zeros( Num, 1 );
potSpot.sl = zeros( Num, 2 );

rg1 = -1:1;   r1 = 2*1+1;
rg2 = -2:2;   r2 = 2*2+1;
rg3 = -3:3;   r3 = 2*3+1;
rg4 = -4:4;   r4 = 2*4+1;
rg5 = -5:5;   r5 = 2*5+1;

i_num(1) = 1;
i_num(2) = 1 / ( r1^3 );
i_num(3) = 1 / ( r1*r2^2 );
i_num(4) = 1 / ( r2*r3^2 );
i_num(5) = 1 / ( r3*r4^2 );
i_num(6) = 1 / ( r4*r5^2 );

i_shl(1) = 1; %Center Pixel 
i_shl(2) = 1/( r1^3 - 1 );          %1 / number pixels in shell 2
i_shl(3) = 1/( r1*r2^2 - r1^3 );    %1 / number pixels in shell 3
i_shl(4) = 1/( r2*r3^2 - r1*r2^2 ); %1 / number pixels in shell 4
i_shl(5) = 1/( r3*r4^2 - r2*r3^2 ); %1 / number pixels in shell 5
i_shl(6) = 1/( r4*r5^2 - r3*r4^2 ); %1 / number pixels in shell 6
Background_median = median( img.layers(:) );
    j=0;
for i=1:Num
    
    x = xyz(i,1);
    y = xyz(i,2);
    z = xyz(i,3);

    pix{1} = img.layers( x, y, z );             % center
    pix{2} = img.layers( x+rg1, y+rg1, z+rg1 ); % cubic Shell 2
	pix{3} = img.layers( x+rg2, y+rg2, z+rg1 ); % Shell 3 
	pix{4} = img.layers( x+rg3, y+rg3, z+rg2 ); % Shell 4 
    pix{5} = img.layers( x+rg4, y+rg4, z+rg3 ); % Shell 5
	pix{6} = img.layers( x+rg5, y+rg5, z+rg4 ); % Shell 6
    
        Sum(1) = pix{1};
        Shl(1) = pix{1};
        Mean(1)= pix{1};
    for k=2:6
        Sum(k) = sum( pix{k}(:) );
        Shl(k) = ( Sum(k) - Sum(k-1) ) * i_shl(k); 
        Mean(k)= Sum(k) * i_num(k);        
    end
    for k=1:5
       Contrast(k) = Mean(k) / Shl(k+1);
       contrast(k) = Shl(k) / Shl(k+1);
    end
            
    potSpot.contrast{i} =  Contrast; 
    

    if sum( Contrast > Thresh  ) >= 3 || ... 
       sum( contrast > 1.04    ) >= 5
 
  
       j=j+1;       
 
       sp.mean(j) =  Mean(4);      
       sp.background(j) = Shl(6);
       Pix = pix{5}(:) - Shl(6);%Background_median;
       sp.intensity(j, 1) = ( Pix' * Pix ) / sum( Pix );
       sp.xyz(j,:) = xyz(i,:);
      %sp.pixs{j} = pix_b;
      %sp.rsq(j) = Rsq;
      %sp.hat_xyz(j,:) = hat_xyz;
       sp.contrast(j) = Shl(6) / Mean(4);
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
fprintf( '4) After Filtering Non Uniqe Spots: \t%d \n', numel(bestUniqe));

if numel( bestUniqe ) > 500
    Intensities = sp.intensity( bestUniqe ); 
    [sorted_intensity, ...
     sorted_intensity_vals ] = sort( Intensities, 'descend' ); 

    bestUniqe = bestUniqe( sorted_intensity_vals(1:400) );
end 
spb.intensity = sp.intensity( bestUniqe );
spb.contrast = sp.contrast( bestUniqe );
%spb.rsq = sp.rsq( bestUniqe );
%spb.xyz = sp.xyz( bestUniqe, : );

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
%[ spb.xyz( spb.km.mRNAs_ind, : ) ];
spb.mat =   spb.intensity( spb.km.mRNAs_ind );


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





















