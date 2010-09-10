Fish View Overview
------------------
generateFishView.py reads results directory output by analyze_experiment_set.m
in modular_image_analysis and creates HTML files useful for viewing the results.

Requirements
------------
Python (tested in 2.6)
Mako template engine (http://www.makotemplates.org/)
Ghostscript (gs) available in the PATH (for PDF to PNG conversion)

Usage
-----
python generateFishView.py path/to/results

