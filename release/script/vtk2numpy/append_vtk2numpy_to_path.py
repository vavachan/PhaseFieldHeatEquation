import os
import os.path

import sys

def append_vtk2numpy_to_path():
    def check_condition_and_do(directory):
        makefile_path = os.path.join(directory, 'Makefile')
        if os.path.exists( makefile_path ):
            with open(makefile_path, 'r') as file:
                for line in file:
                    if line.startswith('OLB_ROOT'):
                        sys.path.append( directory + '/' + line.split(':=')[1].strip() + '/script/vtk2numpy' )
                        return True
        return False

    path = os.getcwd()
    current_dir = os.path.abspath(path)
    root_dir = os.path.abspath(os.sep)
    
    while current_dir != root_dir:
        if check_condition_and_do(current_dir):
            return
        else:
            current_dir = os.path.dirname(current_dir)
    
    # If we reach this point, the path is exhausted without meeting the condition
    raise FileNotFoundError("Path is exhausted without meeting the condition.")
