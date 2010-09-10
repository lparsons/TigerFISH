'''
Created on Sep 9, 2010

@author: lparsons
'''
import os

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
        return experiment_list
    
    def index_file(self):
        return 'index.html'
    

class Experiment(object):
    '''
    classdocs
    '''


    def __init__(self, name, directory, root_directory):
        '''
        Constructor
        '''
        self.name = name
        self.directory = directory
        self.root_directory = root_directory
        self.regions = self.find_regions()
        
    def find_regions(self):
        region_list = []
        for f in os.listdir(self.abs_directory()):
            if os.path.isdir(os.path.join(self.abs_directory(), f)):
                r = Region(f, f)
                region_list.append(r)
        return region_list
    
    def index_file(self):
        return os.path.join(self.directory, 'index.html')
    
    def abs_directory(self):
        return os.path.join(self.root_directory, self.directory)
                
                
class Region(object):
    '''
    classdocs
    '''
    
    def __init__(self, name, directory):
        '''
        Constructor
        '''
        self.name = name
        self.directory = directory
    
    def index_file(self):
        return os.path.join(self.directory, 'index.html')
        