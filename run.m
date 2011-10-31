global params 
output_dir = [ '~/locSpot' filesep 'limit_eth'];   mkdir ( output_dir );
experiment_list = [ output_dir filesep 'experiment_list.txt'];

params.algorithm = '3D'; %{'3D', '2D', '2D_local'};
params.Threshold_Contrast = [1.06  1.1  1.1];
params.Threshold_Intensity = { 0.195, [], []}; %{ 0.22, 0.32, 0.3 };
params.images.path =  '/home/nslavov/data-fish/limit_eth/';
params.images.mask  = 'N*';
params.nmiss = 0;
params.DAPI_Dimensions = 2;
params.FDR_Threshold = 0.01;
params.NULL = [2 1 1]; 
params.Save_Detailed_Results = true;
params.load_results = true; 

% ADJUSTABLE PARAMETERS TO TUNE locSpot TO YOUR IMAGES: 
%params.Region_num = 5;
% If the DAPI contrast in your images differs significantly from the contrast in our images, cell segementation           ls | xargs -I{} mv {}  "E500_yox1_sur7_no{}"
% may not work properly. If that is the case, as can been seen from the interactive output, adjust  params.DAPI_Contrast
% If params.DAPI_Contrast is not defined or empty, locSpot will use the default value of 0.137
params.DAPI_Contrast = 0.137;
	

generate_experiment_list( 2, 5, experiment_list ); %limit_eth=[1:9 11:14] batch_eth=17:35 %batch_glu=25:34;  %341:342 346:353
experiment_set_data = main(experiment_list, output_dir, 'load_results', params.load_results );

system( ['python  ~/locSpot/fish_view/generateFishView.py '  output_dir ] ); 

