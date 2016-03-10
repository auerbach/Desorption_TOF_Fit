# -*- coding: utf-8 -*-

import global_variables_old as gv

def get_book_variable_module_name(module_name):
    module = globals().get(module_name, None)
    book = {}
    if module:
        book = {key: value for key, value in module.__dict__.items() if not (key.startswith('__') or key.startswith('_'))}
    return book


#==============================================================================
# module = globals().get('gv', None)
# 
# dd = module.__dict__
# dd_items = module.__dict__.items()
# 
#==============================================================================
book = get_book_variable_module_name('gv')

for key, value in book.items():
    print("{:<30}{:<100}".format(key, value))
