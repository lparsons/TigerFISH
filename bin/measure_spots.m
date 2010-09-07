function spot_data = measure_spots(spots, image_input, varargin)

ip = inputParser;
ip.FunctionName = 'measure_spots';
ip.addOptional('background', 'local', @ischar);
ip.parse(varargin{:});
global_bg = strcmpi(ip.Results.background, 'global');


%% Read Image File
% Load image file if not already loaded
if (ischar(image_input)) 
    image = load_image(image_input);
else
    image = image_input;
end

layers = size(image.info,1);
bg = zeros(1,layers);
for i=1:layers
    tmp = image.layers(:,:,i);
    bg(i) = median(tmp(:));
end

%% Test each spot
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
% % r5 = 5;  rg5 = -r5:r5;   rr5 = 2*r5+1;
% % r6 = 6;  rg6 = -r6:r6;   rr6 = 2*r6+1;
% 
% i_num1 = 1/( rr1^2 );
% i_num2 = 1/( rr2^2 - rr1^2 ); 
% i_num3 = 1/( rr3^2 - rr2^2 ); 
% i_num4 = 1/( rr4^2 - rr3^2 ); 
% % i_num5 = 1/( rr5^2 - rr4^2 ); 
% % i_num6 = 1/( rr6^2 - rr5^2 ); 

inds = find(spots);
[y x z] = ind2sub(size(image.layers),inds);
yxz = [y x z];

num_spots = size(inds,1);
spot_data_tmp = zeros(num_spots, 4);

for i=1:num_spots
    
    y = yxz(i,1);
    x = yxz(i,2);
    z = yxz(i,3);
    
    % Skip spots too close to the edge of the image
    x_edge = min(x, image.info(1).Width-x);
    y_edge = min(y, image.info(1).Height-y);
    if (x_edge <= spot_radius || y_edge <= spot_radius)
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
% %     pix_5 = image.layers(y+rg5, x+rg5, z);
% %     pix_6 = image.layers(y+rg6, x+rg6, z);
%     
%     rm(1) = sum( pix_1(:) ) * i_num1;
%     rm(2) = (sum( pix_2(:) ) - sum(pix_1(:)) ) * i_num2;
%     rm(3) = (sum( pix_3(:) ) - sum(pix_2(:)) ) * i_num3;
%     rm(4) = (sum( pix_4(:) ) - sum(pix_3(:)) ) * i_num4;
% %     rm(5) = (sum( pix_5(:) ) - sum(pix_4(:)) ) * i_num5;
% %     rm(6) = (sum( pix_6(:) ) - sum(pix_5(:)) ) * i_num6;
    
    if global_bg
        parameters.blacklevel = bg(z);
    else
        parameters.blacklevel = rm(spot_radius);
    end
    %parameters.blacklevel = bg(z);
    parameters.psfwidth = 1.2;
    spot_data_tmp(i,[1,2,4]) = gmask_fit(pix{spot_radius}, parameters);
    spot_data_tmp(i,1) = x-spot_radius + spot_data_tmp(i,1);
    spot_data_tmp(i,2) = y-spot_radius + spot_data_tmp(i,2);
    spot_data_tmp(i,3) = z;
    
    %stats = regstats(rm,1:4,'linear',{'rsquare', 'beta'});
    %disp(stats.rsquare)
    %disp(stats.beta)
    %disp(rm(3)-rm(1))
    %disp(contrast)
    
    %if (x == 33 && y == 424)
    %figure, imshow(mat2gray(pix_4), 'InitialMagnification', 1000);
    %pause
    %end
    
end

%disp(num_spots)
%size(spot_data_tmp);
if ~isempty(spot_data_tmp)
    spot_data = spot_data_tmp(spot_data_tmp(:,1)>0,:);
else
    spot_data = [];
%spot_data = spot_data_tmp;

end