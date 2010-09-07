function spots_cell_idx = map_spots_to_cells(cell_map, spot_xyz, offsets)

[cell_map_labeled,num_cells] = bwlabel(cell_map);
edge_cell_idx = identify_edge_cells(cell_map, []);

num_spots = size(spot_xyz,1);
spots_cell_idx = zeros(num_spots,2);

for i=1:num_spots,
    spot = round(spot_xyz(i,1:3));
    % TODO add offsets to spot_xyz
    spots_cell_idx(i,1) = cell_map_labeled(spot(2),spot(1));
    spots_cell_idx(i,2) = edge_cell_idx(edge_cell_idx(:,1)==spots_cell_idx(i,1),2);
end
end