"""
@author:
lanier4@illinois.edu
mradmstr514226508@gmail.com
"""
import os
import sys
import time
import yaml
from inspect import getmembers, isfunction, getsource, signature, isclass

sys.path.insert(1, '../src');
import z_plane as zp
import graphic_utility as gu;
import itergataters as ig
import numcolorpy as ncp

imported_modules_dict = { 'gu': 'graphic_utility.py'
                        , 'ig': 'itergataters.py'
                        , 'ncp': 'numcolorpy.py'
                        , 'zp': 'z_plane.py'}

conscientious_message = 'USER MISTAKE -- NOT AN ERROR'

def get_source_modules_dict():
    return imported_modules_dict

def display_module(imported_module, show_imported_functions=False):
    """ Usage: display_module(imported_module) 
    Args:
        imported_module:         an imported python module
        show_imported_functions: default is False - ignore imported functions
    """
    ignore_functions_list = ['getmembers', 'isfunction', 'getsource', 'signature', 'isclass']
    
    #                                                            display classes in the module
    class_list = [o for o in getmembers(imported_module) if isclass(o[1])]
    if len(class_list) >= 1:
        print('\n\t\tClasses in module %s\n'%(imported_module.__name__))
        
    for class_list_number in range(len(class_list)):
        selected_class = class_list[class_list_number][1]
        class_members = getmembers(selected_class)
        class_methods = [o for o in class_members if isfunction(o[1])]
        print('\nclass %s\n'%(class_list[class_list_number][0]))
        try:
            selected_class_source = getsource(selected_class)
            docs_string = selected_class_source.split('"""')[1]
            if docs_string is not None and len(docs_string) > 0:
                print(docs_string, '\n')
        except:
            pass
        
        display_methods_functions_list(class_methods, string_prefix='\t')
                
    #                                                            display functions in the module
    functions_list = [o for o in getmembers(imported_module) if isfunction(o[1])]
    if len(functions_list) >= 1:
        print('\n\t\tfunctions in module %s\n'%(imported_module.__name__))

    display_methods_functions_list(functions_list, 
                                   string_prefix='', 
                                   ignore_functions_list=ignore_functions_list,
                                   show_imported_functions=True)
                
def display_methods_functions_list(meth_func_list, 
                                   string_prefix='',
                                   ignore_functions_list=[],
                                   show_imported_functions=False):
    source_list = [getsource(o[1]) for o in meth_func_list if isfunction(o[1])]
    if len(source_list) != len(meth_func_list): return  #       This should not be possible
    
    for list_number in range(len(meth_func_list)):
        function_tuple = meth_func_list[list_number]    
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

            print(string_prefix, 'def %s%s'%(function_tuple[0], signature(function_tuple[1])))
            if docs_string is None:
                print('doc_missing\n')
            else:
                print(docs_string,'\n')