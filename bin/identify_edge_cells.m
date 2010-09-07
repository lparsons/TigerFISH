function cell_idx = identify_edge_cells(cell_map, offsets)

[cell_map_labeled,num_cells] = bwlabel(cell_map);

nonbordercells = imclearborder(cell_map_labeled);
nonbordercells = unique(nonbordercells);
% nonbordercells = nonbordercells(nonbordercells~=0);

% cell_idx = true(num_cells+1,1);
% cell_idx(nonbordercells) = false;

cell_idx = ones(numel(unique(cell_map_labeled)),2);
cell_idx(:,1) = unique(cell_map_labeled);
cell_idx(nonbordercells+1,2) = 0;
cell_idx(cell_idx(:,1)==0,2) = 2;

% % Remove cells touching the edge of the image
% % Create padding based on image registration
% paddedCellMap = ones(size(allCellMap));
% if (registration_output(3) > 0)
%     padding(1) = ceil(registration_output(3));
%     padding(3) = 0;
% else
%     padding(1) = 0;
%     padding(3) = -floor(registration_output(3));
% end
% if (registration_output(4) > 0)
%     padding(2) = ceil(registration_output(4));
%     padding(4) = 0;
% else
%     padding(2) = 0;
%     padding(4) = -floor(registration_output(4));
% end
% %trimRegion = [padding(1):size(allCellMap, 1)-padding(3), padding(2):size(allCellMap, 2)-padding(4)];
% trimedCellMap = allCellMap(1+padding(1):size(allCellMap, 1)-padding(3), 1+padding(2):size(allCellMap, 2)-padding(4));
% paddedCellMap(1+padding(1):size(allCellMap, 1)-padding(3), 1+padding(2):size(allCellMap, 2)-padding(4)) = trimedCellMap;
% 
% % Remove cells connected to border plus padding
% nonbordercells = imclearborder(paddedCellMap);
% bordercells = allCellMap & ~nonbordercells;

end