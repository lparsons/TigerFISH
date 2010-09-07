% gmask_fit calculates the center of a psf using an iterative gaussian mask (Thompson et al.)
%
% Based on gmask_fit IDL program by Dan Larson - 8/16/00
% MatLab implementation - Lance Parsons, 4/15/2009


function results = gmask_fit(pic, parameters)

x_dim = size(pic,2);
y_dim = size(pic,1);
x0=0;                               % center of mass coordinates
y0=0;
F=1/(sqrt(2)*parameters.psfwidth);  % This factor shows up repeatedly in GMASK
gauss_mask = zeros(size(pic));      % Gaussian mask for centroid fitting
error = 0;                          % RMS distance between the 'true' spot and the centroid spot
results = zeros(3,1);                 % this array holds the returned results of the function:
%                                   % x0, y0, number of photons in spot

% offset correction
blksubtract=(pic-parameters.blacklevel);
% Set floor to zero (no negative values)
image=(blksubtract>0).*blksubtract;

% boundary condition.  border is set to zero
image(1, :)= 0;
image(x_dim, :)= 0;
image(:, 1)= 0;
image(:, y_dim) = 0;

% easy localisation by finding the centroid of the image
spot_pixels = ones(size(image));
spot_pixels(image==0) = 0;
stats = regionprops(spot_pixels, image, 'Centroid');
if (size(stats,1)~=0)
    center=stats.Centroid;
    x0=center(1);
    y0=center(2);
else
    results(1)=floor(x_dim/2);
    results(2)=floor(y_dim/2);
    results(3)=0.0;
    return
end



%% iterative centroid calculation with gaussian mask
h = 1.0e-3; %tolerance
diff_x=0.0;
diff_y=0.0;
repeat_index=0;

arr=reshape(0:(x_dim*y_dim)-1, x_dim,y_dim);
xarr=mod(arr,x_dim);
yarr=arr/x_dim;


while (((abs(diff_x) >= h) && (abs(diff_y) >= h)) && (repeat_index < 300)) || repeat_index == 0
    x0=x0 + diff_x/2.0;
    y0=y0 + diff_y/2.0;
    a=F*(yarr - 0.5 - y0);
    b=F*(yarr + 0.5 - y0);
    c=F*(xarr - 0.5 - x0);
    d=F*(xarr + 0.5 - x0);
    gauss_mask =0.25.*(erf(a)-erf(b)).*(erf(c)-erf(d));
    total=sum(image.*gauss_mask);
    trial_x0=sum(xarr.*image.*gauss_mask);
    trial_y0=sum(yarr.*image.*gauss_mask);
    %print ([x0, trial_x0, total]);
    %print trial_x0/total;
    diff_x = trial_x0/total - x0;
    diff_y = trial_y0/total - y0;
    repeat_index=repeat_index+1;
end

if (repeat_index > 300)
    sprintf('GMASK ITERATION MAXED OUT (number of iterations=%d)\n', repeat_index)
    results(3)=0.0;
    return
end

%sprintf('gmask_fit convergence satisfied (number of iterations=%d)\n', repeat_index)

% photon number calc
total = sum(gauss_mask.*gauss_mask);
N = sum(image.*gauss_mask);
photon_number = N/total;

if x0==0
    figure, imshow(mat2gray(pic),'InitialMaginification', 500)
    pause
end

results(1)=x0;
results(2)=y0;
results(3)=photon_number;
return

end