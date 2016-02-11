# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 15:50:15 2016

@author: dja
"""
import pyreadline as readline


#==============================================================================
# def rlinput(prompt, prefill=''):
#    readline.set_startup_hook(lambda: readline.insert_text(prefill))
#    try:
#       return input(prompt)
#    finally:
#       readline.set_startup_hook()
#       
# if __name__ == '__main__':
#     name = rlinput('file_name: ', prefill='file1.txt')
#     print()
#     print(name)
#==============================================================================


#==============================================================================
# def startup_hook():
#     readline.insert_text('from startup_hook')
# 
# def pre_input_hook():
#     readline.insert_text(' from pre_input_hook')
#     readline.redisplay()
# 
# readline.set_startup_hook(startup_hook)
# readline.set_pre_input_hook(pre_input_hook)
# readline.parse_and_bind('tab: complete')
# 
# while True:
#     line = input('Prompt ("stop" to quit): ')
#     if line == 'stop':
#         break
#     print('ENTERED: "%s"' % line)
#==============================================================================
    
    
    
import win32console

_stdin = win32console.GetStdHandle(win32console.STD_INPUT_HANDLE)

def input_def(prompt, default='junk'):
    keys = []
    for c in str(default):
        evt = win32console.PyINPUT_RECORDType(win32console.KEY_EVENT)
        evt.Char = c
        evt.RepeatCount = 1
        evt.KeyDown = True
        keys.append(evt)

    _stdin.WriteConsoleInput(keys)
    return input(prompt)

if __name__ == '__main__':
    name = input_def('Folder name: ')
    print()
    print(name)