Quick usage guide
-----------------

1) Prepare experiment set data structure
	a) User wrapper.m to parse 'standard' directory structures produced by Ryan
		experiment_set = wrapper(Path, Experiment_Numbers, [filemask=*], [output_dir=output])
	b) Use parse_experiment_dir.m to parse a 'flat' directory (as in 96 well test dir)
		experiment_set = parse_experiment_dir(exp_dir, [filemask='*'], [region_marker='Position']);

2) Analyze experiment set (currently only Lance's method implemented)
	a) experiment_set_data = analyze_experiment_set(experiment_set, 'output')

3) Run generateFishView.py to generate HTML viewing pages
    a) python generateFishView.py path/to/results -n 'Experiment Set Name'