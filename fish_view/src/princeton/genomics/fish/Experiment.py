'''
Created on Sep 9, 2010

@author: lparsons
'''

import glob
import os
import re

class Experiment_Set(object):
    '''
    classdocs
    '''
    
    def __init__(self, name, directory):
        '''
        Constructor
        '''
        self.name = name
        self.directory = directory
        self.experiments = self.find_experiments()
        
    def find_experiments(self):
        experiment_list = []
        for f in os.listdir(self.directory):
            if os.path.isdir(os.path.join(self.directory, f)):
                e = Experiment(f, f, self.directory)
                experiment_list.append(e)
            experiment_list.sort(key=lambda x: x.number)
        return experiment_list
    
    def index_file(self):
        return 'index.html'
    

class Experiment(object):
    '''
    Create with name and directory
    Number is set to name.split('_')[1] if that is a float, else it is name 
    '''

    def __init__(self, name, directory, root_directory):
        '''
        Constructor
        '''
        self.name = name
        self.directory = directory
        self.root_directory = root_directory
        self.number = self.name
        name_match = re.match(
            "(?P<experiment>[^_]+)_(?P<gene1>[^_]+)_(?P<gene2>[^_]+)_(?P<gene3>[^_]+)(_(?P<notes>[^_]+))*",
            self.name)
        if name_match:
            num_match = re.match("(?P<text>[a-zA-Z]+)*(?P<number>[\d\.]+)",
                         name_match.group('experiment'))
            if num_match:
                try:
                    self.number = float(num_match.group('number'))
                except:
                    pass
        self.regions = self.find_regions()
        
    def find_regions(self):
        region_list = []
        for f in os.listdir(self.abs_directory()):
            if os.path.isdir(os.path.join(self.abs_directory(), f)):
                r = Region(f, f)
                region_list.append(r)
            region_list.sort(key=lambda x: x.number)
        return region_list
    
    def index_file(self):
        return os.path.join(self.directory, 'index.html')
    
    def abs_directory(self):
        return os.path.join(self.root_directory, self.directory)
    
    def histogram(self, dye, format='png', force=False):
        filename = ''
        if format == 'pdf':
            filename = '%s_spot_intensity_histogram.pdf' % dye
        elif format == 'png':
            filename = '%s_spot_intensity_histogram.png' % dye
            if (not os.path.exists(os.path.join(self.abs_directory(), filename))) or force:
                cmd = "gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -dTextAlphaBits=4 -sOutputFile='%s' '%s'" % (os.path.join(self.abs_directory(), filename), os.path.join(self.abs_directory(), self.histogram(dye, 'pdf')))
                print cmd
                os.system(cmd)
        return filename
    
    def joint_distributions(self, type='prob', force=False):
        filenames = []
        file_paths = glob.glob(os.path.join(self.abs_directory(), 'joint_dist_%s_*.pdf' % type))
        for f in file_paths:
            f_png = f.replace('.pdf', '.png')
            if (not os.path.exists(f_png)) or force:
                cmd = "gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -dTextAlphaBits=4 -sOutputFile='%s' '%s'" % (f_png, f)
                print cmd
                os.system(cmd)
            filenames.append((os.path.basename(f_png), os.path.basename(f)))
        return filenames
            
    def dna_content(self, force=False):
        f_pdf = 'DNA_content.pdf'
        f_png = 'DNA_content.png'
        if (not os.path.exists(os.path.join(self.abs_directory(), f_png))) or force:
            cmd = "gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -dTextAlphaBits=4 -sOutputFile='%s' '%s'" % (os.path.join(self.abs_directory(), f_png), os.path.join(self.abs_directory(), f_pdf))
            print cmd
            os.system(cmd)
        return (f_png, f_pdf)
                
class Region(object):
    '''
    Create with name and directory
    Number is set to name.split('_')[1] if that is a float, else it is name 
    '''
    
    def __init__(self, name, directory):
        '''
        Constructor
        '''
        self.name = name
        self.directory = directory
        self.number = name
        try:
            number = name.split('_')[1]
            self.number = float(number)
        except:
            pass

    def index_file(self):
        return os.path.join(self.directory, 'index.html')
        
