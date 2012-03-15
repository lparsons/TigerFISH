function [ new_spot_data ] = merge_spots( varargin )
%MERGE_SPOTS Summary of this function goes here
%   Detailed explanation goes here
ip = inputParser;
ip.FunctionName = 'main';
ip.addRequired('spot_data',@isnumeric);
ip.addRequired('duplicateThreshold', @isnumeric);
ip.parse(varargin{:});
spot_data = ip.Results.spot_data;
duplicateThreshold = ip.Results.duplicateThreshold;

if size(spot_data,1) <= 1
    new_spot_data = spot_data;
else
    
    % Sort based on intensity
    [b, ix] = sort(spot_data(:,4));
    % Get rank for each spot
    ix2rank = [ix, (1:length(ix))'];
    % Assign rank to each spot
    temp_spot_data = spot_data;
    temp_spot_data(ix(:,1),6) = ix2rank(:,2);
    % Consolidate, choosing the row with the highest intensity
    % when there are duplicates within tolerance
    % Note: Z is multiplied by .25 so duplicates in nearby layers are
    % removed, must be at least a few layers apart
    [xc,yc] = consolidator13([spot_data(:,1:2), spot_data(:,3)*.25],temp_spot_data(:,6),'max',duplicateThreshold);
    % Select only the consolidated spots
    unique_spots = ismember(temp_spot_data(:,6),yc);
    new_spot_data = spot_data(unique_spots,:);
end

end

