#!/usr/bin/env python

'''
Created on Sep 9, 2010

@author: lparsons
'''
import os
import sys
sys.path.append(os.path.join(sys.path[0], "src"))

from mako import exceptions
from mako.lookup import TemplateLookup
from mako.template import Template
from optparse import OptionParser
from princeton.genomics.fish.Experiment import Experiment_Set




def main():
    usage = "usage: %prog [options] results_directory"    
    
    parser = OptionParser(usage=usage)
    parser.add_option("-n", "--name", dest="set_name", type="string", help="Name of the Experiment Set")
    parser.add_option("-f", "--force", dest="force", action="store_true",
                      default=True, help="Force regeneration of png files from "
                      "pdf files (default: %default)")
    (opts, args) = parser.parse_args()
    
    if len(args) < 1:
        parser.print_help()
        sys.exit()
        
    results_directory = args[0]
    
    if (opts.set_name):
        set_name = opts.set_name
    else:
        set_name = os.path.basename(results_directory)
        
    experiment_list = Experiment_Set(set_name, results_directory)

    # Experiment index
    mylookup = TemplateLookup(directories=[os.path.join(sys.path[0], 'templates')], module_directory=os.path.join(sys.path[0], 'tmp/mako_modules'))
    index_file = open(os.path.join(results_directory, 'index.html'), 'w')
    try:
        experiment_set_template = mylookup.get_template('experiment_set.txt')
        index_file.write(experiment_set_template.render(experiment_set=experiment_list))
    except:
        print exceptions.text_error_template().render()
    index_file.closed
    
    # Experiment View
    for e in experiment_list.experiments:
        experiment_index_file = open(os.path.join(e.root_directory, e.index_file()), 'w')
        try:
            experiment_template = mylookup.get_template('experiment.txt')
            experiment_index_file.write(experiment_template.render(experiment=e, force=opts.force))
        except:
            print exceptions.text_error_template().render()
        
        # Region View
        for r in e.regions:
            region_index_file = open(os.path.join(e.abs_directory(), r.directory, 'index.html'), 'w')
            try:
                region_template = mylookup.get_template('region.txt')
                region_index_file.write(region_template.render(experiment=e, region=r, force=opts.force))
            except:
                print exceptions.text_error_template().render()

    


        
if __name__ == '__main__':
    main()
    sys.exit()
    
