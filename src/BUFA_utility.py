"""
@author:
lanier4@illinois.edu
mradmstr514226508@gmail.com
"""

import os
import time
import yaml
from inspect import getmembers, isfunction, getsource, signature, isclass

conscientious_message = 'USER MISTAKE -- NOT AN ERROR'

def display_module(imported_module, show_imported_functions=False):
    """ Usage: display_module_functions(any_module) 
    Args:
        imported_module:         an imported python module
        show_imported_functions: default is False - ignore imported functions
        
    """
    ignore_functions_list = ['getmembers', 'isfunction', 'getsource', 'signature']

    functions_list = [o for o in getmembers(imported_module) if isfunction(o[1])]
    source_list = [getsource(o[1]) for o in getmembers(imported_module) if isfunction(o[1])]    
    
    if len(source_list) != len(functions_list): return  #       This should not be possible
    
    for list_number in range(len(functions_list)):
        function_tuple = functions_list[list_number]    
        if function_tuple[0] in ignore_functions_list:
            if show_imported_functions == True:
                print('using: %s%s'%(function_tuple[0], signature(function_tuple[1])))
        else:
            source_str = source_list[list_number]
            docs_string = None
            try:
                docs_string = source_str.split('"""')[1]
            except:
                pass

            print('def %s%s'%(function_tuple[0], signature(function_tuple[1])))
            if docs_string is None:
                print('doc_missing\n')
            else:
                print(docs_string,'\n')
                
                
    class_list = [o for o in getmembers(imported_module) if isclass(o[1])]
    class_source_list = [getsource(o[1]) for o in getmembers(imported_module) if isclass(o[1])]
    for list_number in range(len(class_list)):
        class_tuple = class_list[list_number]    
        source_str = class_source_list[list_number]
        docs_string = None
        try:
            docs_string = source_str.split('"""')[1]
        except:
            pass

        print('def %s'%(class_tuple[0]))
        if docs_string is None:
            print('doc_missing\n')
        else:
            print(docs_string,'\n')