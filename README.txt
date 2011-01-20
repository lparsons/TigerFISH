Usage Guide
-----------

1) Prepare an tab delimted file listing the experiments an inputfiles

    Coumns: Experiment, Region, Cy3_label, Cy3_file, Cy3.5_label, Cy3.5_file, Cy5_label, Cy5_file, DAPI_label, DAPI_file

    a) This can be generated using the wrapper function for files in a directory structure like:
        EXPERIMENT_SET/
            EXPERIMENT/
                1/
                    EXPERIMENT_NUMBER-EXPERIMENT_NAME - Position #_CDYE.tiff
                        DYE should be one of CY3, CY3.5, CY5, DAPI
                        e.g. NS1-POL30_SUR4_OM45 - Position 1_CCY3.tiff

        genereate_experiment_list(path, 
                                  [experiment_mumbers=1:1000],
                                  [output_filename='experiment_list.txt'], 
                                  [filemask=*])

2) Analyze experiment set
    main( experiment_list_file, output_dir, [algorithm='3D'], [load_results=False] )

3) Run generateFishView.py to generate HTML viewing pages
    a) python generateFishView.py path/to/results -n 'Experiment Set Name'




OLD METHOD (Deprecated)
-----------------------
1) Prepare experiment set data structure
	a) User wrapper.m to parse 'standard' directory structures produced by Ryan
		experiment_set = wrapper(Path, Experiment_Numbers, [filemask=*], [output_dir=output])
	b) Use parse_experiment_dir.m to parse a 'flat' directory (as in 96 well test dir)
		experiment_set = parse_experiment_dir(exp_dir, [filemask='*'], [region_marker='Position']);

2) Analyze experiment set (currently only Lance's method implemented)
	a) experiment_set_data = analyze_experiment_set(experiment_set, 'output')

3) Run generateFishView.py to generate HTML viewing pages
    a) python generateFishView.py path/to/results -n 'Experiment Set Name'