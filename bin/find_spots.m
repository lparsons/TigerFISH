function spots = find_spots(image_input, contrast_threshold)


%% Read Image File
% Load image file if not already loaded
if (ischar(image_input))
    image = load_image(image_input);
else
    image = image_input;
end

%if nargin <2, contrast_threshold = 1.05; end 
%if nargin <3, Rsq_Treshold = 0.4; end 


%% Identify potential spots
% if ~isempty(regexpi(image.info(1).ImageDescription, 'EMGain', 'match'))
%     % New image, no denoising
image = tophat_image(image);
%J = image.tophat_max;
% else
%     % Old image, denoiseing required
%     tmp = wiener2(mat2gray(image.max));
%     se = strel('square',5);
%     J = imtophat(tmp, se);
% end
J = mat2gray(image.tophat_max);
thresh = median(J(:));
mask = im2bw(J, thresh*2);
mask = bwareaopen(mask, 2);
%mask = imextendedmax(J .* mask, .1);
mask = imregionalmax(J .* mask);

% Debug
%ovr = imoverlay(mat2gray(image.max), mask, [1 0 0]);
%figure, imshow(ovr);


%% ????
%[spots, num_spots] = bwlabel(mask);
%bb = regionprops(spots, image.max, 'BoundingBox', 'PixelList');


%region = image.layers(bb(1).PixelList(:,1), bb(1).PixelList(:,2),:);
%[zmax.vals  zmax.inds] = max(region, [], 3);


%% Find max layer for each pixel
% [zmax.vals  zmax.inds] = max(image.layers, [], 3);
% [y x] = find(mask);
% z = zmax.inds(mask);
%
% spots = false(size(image.layers));
% spots(sub2ind(size(image.layers),y,x,z)) = true;


% xyz = [x y z];
% num_spots = size(xyz, 1);
%
%% Test each spot
spots = false(size(image.layers));
[regions, num_spots] = bwlabel(mask);
spot_props = regionprops(regions, image.max, 'BoundingBox', 'PixelList', 'Centroid');

% Setup some variables to get pixels around spot center
spot_radius = 4;
rg = cell(1,spot_radius); rr = zeros(1,spot_radius); i_num = zeros(1, spot_radius);
for r=1:spot_radius
    rg{r} = -r:r;
    rr(r) = 2*r+1;
    i_num(r) = 1/( rr(r)^2 );
    if (r>1)
        i_num(r) = 1/(rr(r)^2-rr(r-1)^2);
    end
end

% r1 = 1;  rg1 = -r1:r1;   rr1 = 2*r1+1;
% r2 = 2;  rg2 = -r2:r2;   rr2 = 2*r2+1;
% r3 = 3;  rg3 = -r3:r3;   rr3 = 2*r3+1;
% r4 = 4;  rg4 = -r4:r4;   rr4 = 2*r4+1;
%
% i_num1 = 1/( rr1^2 );
% i_num2 = 1/( rr2^2 - rr1^2 );
% i_num3 = 1/( rr3^2 - rr2^2 );
% i_num4 = 1/( rr4^2 - rr3^2 );

for i=1:num_spots
    
    bb = spot_props(i);
    [xy] = floor(bb.Centroid);
    x_cent = xy(1);
    y_cent = xy(2);
    %
    %     x = xyz(i,1);
    %     y = xyz(i,2);
    %     z = xyz(i,3);
    %
    % Skip if too close to edge of image
    x_edge = min(x_cent, image.info(1).Width-x_cent);
    y_edge = min(y_cent, image.info(1).Height-y_cent);
    if (x_edge <= spot_radius || y_edge <= spot_radius)
        continue
    end
    
    % Find max layer
    %     max_pix_3 = image.max(y_cent+rg{3}, x_cent+rg{3});
    %     [mval mind] = max(max_pix_3(:));
    %     [my mx] = ind2sub(size(max_pix_3), mind);
    max_pix_1 = image.max(y_cent+rg{1}, x_cent+rg{1});
    [mval mind] = max(max_pix_1(:));
    [my mx] = ind2sub(size(max_pix_1), mind);
    
    x = x_cent-1 + mx -1;
    y = y_cent-1 + my -1;
    [mzval z] = max(image.layers(y,x,:));
    
    % Skip if too close to edge of image
    x_edge = min(x, image.info(1).Width-x);
    y_edge = min(y, image.info(1).Height-y);
	z_edge = min(z, size(image.layers,3)-z);
    if (x_edge <= spot_radius || y_edge <= spot_radius ) %|| z_edge <= spot_radius)
        continue
    end
    
    pix = cell(1,spot_radius); rm = zeros(1, spot_radius);
    for r=1:spot_radius
        pix{r} = image.layers(y+rg{r}, x+rg{r}, z);
        if r==1
            rm(r) = sum( pix{r}(:) ) * i_num(r);
        else
            rm(r) = (sum( pix{r}(:) ) - sum(pix{r-1}(:)) ) * i_num(r);
        end
    end
    
    %     pix_1 = image.layers(y+rg1, x+rg1, z);
    %     pix_2 = image.layers(y+rg2, x+rg2, z);
    %     pix_3 = image.layers(y+rg3, x+rg3, z);
    %     pix_4 = image.layers(y+rg4, x+rg4, z);
    %
    %     rm(1) = sum( pix_1(:) ) * i_num1;
    %     rm(2) = (sum( pix_2(:) ) - sum(pix_1(:)) ) * i_num2;
    %     rm(3) = (sum( pix_3(:) ) - sum(pix_2(:)) ) * i_num3;
    %     rm(4) = (sum( pix_4(:) ) - sum(pix_3(:)) ) * i_num4;
	
	
    
    valid_spot = true;
    contrast = rm(1)/rm(spot_radius);
    if (contrast > contrast_threshold)
        for r=1:(spot_radius-1)
            if rm(r) < rm(r+1)
                valid_spot = false;
            end
        end
    else
        valid_spot = false;
    end
    
    if valid_spot
        spots(y,x,z) = true;
		
%  		[sl Mean Rsq] = spot_regress_fun(  image.layers(y+rg{3}, x+rg{3}, z+rg{3}), 7 ); %Rsq, pause (0.5)
%  		if Rsq > Rsq_Treshold
%  		  spots(y,x,z) = true;
%  		end
    end
	
	
    
    %     if (rm(1) > rm(2) && rm(2) > rm(3) && contrast > 1.1)
    %         spots(y,x,z) = true;
    %     end
    
    
    %     %stats = regstats(rm,1:4,'linear',{'rsquare', 'beta'});
    %     %disp(stats.rsquare)
    %     %disp(stats.beta)
    %     %disp(rm(3)-rm(1))
    %     %disp(contrast)
    %
    %     %if (x == 33 && y == 424)
    %     %figure, imshow(mat2gray(pix_4), 'InitialMagnification', 1000);
    %     %pause
    %     %end
    %
end