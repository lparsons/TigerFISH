function [ image ] = tophat_image( varargin )
%tophat_image Summary of this function goes here
%   Detailed explanation goes here

[image, se] = parse_inputs(varargin{:});


% Read image layers
image.tophat_layers = zeros(image.info(1).Height, image.info(1).Width, size(image.info, 1));
for i=1:1:size( image.info, 1 )
    image.tophat_layers(:,:,i) = imtophat(image.layers(:,:,i), se);
end

% Compute max projection and store mid-slice
image.tophat_max = max(image.tophat_layers, [], 3);
image.tophat_mid = image.tophat_layers(:,:,floor(size(image.info, 1)/2));


%%%
%%% parse_inputs
%%%
    function [image, se] = parse_inputs(varargin)
        
        iptchecknargin(1,2,nargin,mfilename)
        
        image = varargin{1};
        iptcheckinput(image,{'struct'},{'nonempty'},mfilename,'image',1);
        
        if (nargin >= 2)
            se = varargin{2};
        else
            se = strel('square',5);
        end
        iptcheckinput(se,{'strel'},{'nonempty'},mfilename,'se',2);
    end

end