#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module to parse the control file for FitDesorbtionTOF

"""

__author__ = 'Daniel Auerbach'  
__email__  = 'daniel@djauerbach.com'

from glob import glob
from Parameters2 import Parameters2
import TOF_fit_global

class Fit_control(object):
    def __init__(self):
        # names of parameters constants and lists
        self.parms_list = [
                           'e0',
                           'w',
                           'tcutc', 'tcutw',
                           'ecutm', 'ecuts',
                           'temp',
                           'ffr',
                           'ion_tof',
                           'y_scale',
                           'baseline']

        self.vars_list = ['comment_xlsx',
                          'grid_type',
                          'l_laser',
                          'points_detector',
                          'points_source',
                          'r_aperture',
                          'r_final',
                          'r_source',
                          'r_diff_wall',
                          'theta_step',
                          'z_aperture',
                          'z_diff_wall',
                          'z_final',
                          'z_laser',
                          'z_refef',
                          'z_source'
                          ]

        self.lists_list = ['averaging_type',
                           'cutoff_type',
                           'function',
                           'n_delt',
                           'threshold'
                           ]

        self.tuples_list = ['fit_range',
                            'baseline_range']

        self.errors = []
        self.global_parms=[]
        self.parms = None


    #==========================================================================
    # print list items
    #==========================================================================
    def print_list_items(self):
        for item in self.lists_list:
            # print('{:16s} = {}'.format(item + 's', eval('glbl_old.' + item + 's')))
            print('{:16s} = {}'.format(item + 's', eval('glbl.' + item + 's')))


    #==========================================================================
    # print tuple items
    #==========================================================================
    def print_tuple_items(self):
        for item in self.tuples_list:
            print('{:16s} = {}'.format(item + 's', eval('glbl.' + item + 's')))


    #==============================================================================
    # open file for reading; abort if file not found
    #==============================================================================
    def openForReading(self, filename):
        try:
            file=open(filename, 'r')
            return file
        except FileNotFoundError:
            print('\n ***')
            print(' *** File ', filename, ' not found')
            print(' *** Ending program')
            print(' *** \n')
            raise SystemExit


    #==============================================================================
    # break input line into tokens
    #=========================================================================
#==============================================================================
#     def getTokens(self, line):
# 
#         # split the line on '=,'
#         cmd = line.split('=')[0].strip()
#         
#         try:
#             arg = line.split('=')[1].strip()
#             tokens = arg.split(',')
#         except:
#             arg = None
#             tokens = None
#         
#         if tokens:
#             for j in range(len(tokens)):
#                 tokens[j] = tokens[j].strip()
#                 # replace nule tokens (i.e. '') with None
#                 if tokens[j] == '': tokens[j] = None
#                 # remove imbedded comments
#                 try:
#                     ii = tokens[j].index('#')
#                     # if so, remove rest of token and strip white space
#                     tokens[j] = tokens[j][0:ii].strip()
#                     # ignore rest of tokens
#                     for k in range(j + 1, len(tokens)):
#                         tokens.pop(-1)
#                     break
#                 except:
#                     pass
# 
#         return cmd, arg, tokens
#==============================================================================
    
    def tokenize(self, line, sep, count):
        toks = line.split(sep, count)
        for j in range(len(toks)):
            toks[j] = toks[j].strip()
            
            # replace nule tokens (i.e. '') with None
            if toks[j] == '': 
                toks[j] = None
            
            # remove imbedded comments
            # try to get index of # symbol in toks[j]
            # if found remove rest of token and strip white space
            try:
                ii = toks[j].index('#')     
                toks[j] = toks[j][0:ii].strip()
                # ignore rest of tokens
                for k in range(j + 1, len(toks)):
                    toks.pop(-1)
                break
            except:
                pass
        return toks

        
        
    def getTokens(self, line):

        # split the line on '=,'
        toks = self.tokenize(line, '=', 1)
        cmd = toks[0].lower()
        
        # if there is more than 1 token break apart further
        # sometimes there will only be one as in 'End Section'
        if len(toks)>1:
            arg = toks[1].lower()
            # split the args on ',' to get tokens
            tokens = self.tokenize(arg, ',', -1)
        else:
            arg = None
            tokens = None

        return cmd, arg, tokens


    #==============================================================================
    # check if cmd is a Parameter item
    #==============================================================================
    # def isParameter(self, tokens):
    #     return any(tokens[0] in s for s in self.parms_list)
    def isParameter(self, cmd):
        return any(cmd in s for s in self.parms_list)


    #==============================================================================
    # check if cmd is an Global Variable item
    #==============================================================================
    def isVar(self, cmd):
        return any([cmd == s for s in self.vars_list])


    #==============================================================================
    # check if cmd is a List item
    #==============================================================================
    def isListItem(self, cmd):
        matching = [cmd == s for s in self.lists_list]
        return any(matching)

    #==============================================================================
    # check if cmd is a Tuple item
    #==============================================================================
    def isTupleItem(self, cmd):
        matching = [cmd == s for s in self.tuples_list]
        return any(matching)


    #==============================================================================
    # check item is a number
    #==============================================================================
    def isNumber(self, s):
        try:
            float(s)
            return True
        except:
            return False


    #==============================================================================
    # add Parameter
    #==============================================================================
    def addParameter(self, runNumber, cmd, arg, tokens, line, lineNumber):

        # global errors,  global_parms, parms

        # check if number of tokesn is correct
        #   there should be name, and 2, 3 or 4 parameters
        #   i.e. 3 to 5 tokens on the line
        # if len(tokens) < 3 or len(tokens) > 5:
        #     print('\n *** wrong number of items on line ', lineNumber)
        #     print(' *** items found =', len(tokens),
        #           '. Should be 3, 4, or 5')
        #     print(' *** line =', line)
        #     print(' *** items =', tokens)
        #     self.errors.append(('wrong number of items', 'line ' + str(lineNumber)))
        #     return

        if len(tokens) < 2 or len(tokens) > 4:
            print('\n *** wrong number of values for Parameter Item on line ', lineNumber)
            print(' *** number found =', len(tokens),
                  '. Should be 2, 3, or 4')
            print(' *** line =', line)
            print(' *** items =', tokens)
            self.errors.append(('wrong number of items', 'line ' + str(lineNumber)))
            return


        # set defaults
        vary1 = False
        glbl1 = False
        min1 = None
        max1 = None

        # # set vary and global
        # if tokens[2] != '0':
        #         vary1 = True
        # if tokens[2] == '2':
        #     glbl1 = True
        #     self.global_parms.append(tokens[0])

        # set vary and global
        if tokens[1] != '0':
                vary1 = True
        if tokens[1] == '2':
            glbl1 = True
            self.global_parms.append(cmd)

        # set bounds
        # try:
        #     min1=float(tokens[3])
        # except:
        #     pass
        # try:
        #     max1 = float(tokens[4])
        # except:
        #     pass
        try:
            min1=float(tokens[2])
        except:
            pass
        try:
            max1 = float(tokens[3])
        except:
            pass

        self.parms.add(cmd + '_' + str(runNumber),
                  value = float(tokens[0]),
                  vary = vary1,
                  glbl = glbl1,
                  min = min1,
                  max = max1
                  )
        # print('addParameter runNumber, tokens =', runNumber, tokens)


    #==============================================================================
    # add variable
    #==============================================================================
    def addVar(self, cmd, arg, tokens, line, lineNumber, glbl):
        if len(tokens) != 1:
            error = 'wrong number of args for item on line ' + str(lineNumber)
            print('\n ***' + error, lineNumber)
            print(' *** items found =', len(tokens),
                  '. Should be 1')
            print(' *** line =', line)
            print(' *** items =', tokens)
            self.errors.append(error)
            return

        exec('glbl.' + cmd + '=' + tokens[0])


    #==============================================================================
    # add list item
    #==============================================================================
    def addListItem(self, cmd, arg, tokens, line, lineNumber, glbl):

        # global errors,  global_parms, parms

        if len(tokens) != 1:
            error = 'wrong number of arguments items on line ' + str(lineNumber)
            self.errors.append(error)

            print('\n ***' + error)
            print(' *** items found =', len(tokens),
                  '. Should be 1')
            print(' *** line =', line)
            print(' *** items =', tokens)


        # item is a string, need to add quotes to exec, e.g. ERF become 'ERF'
        if self.isNumber(tokens[0]):
            exec('glbl.' + cmd + '=' + tokens[0])
        else:
            exec('glbl.' + cmd + '=' + "tokens[0]")

    #==============================================================================
    # add tuple item
    #==============================================================================
    def addTupleItem(self, cmd, arg, tokens, line, lineNumber, glbl):
        # put tokens[1:] into single string
#        toks = tokens[0]
#        for tok in tokens[1:]:
#            if not tok:
#                tok =''
#            toks += ', ' + tok
#
#        # item is a string, need to add quotes to exec, e.g. ERF become 'ERF'
#        if self.isNumber(tokens[0]):
#            exec('glbl.' + cmd + '=' + toks)
#        else:
#            exec('glbl.' + cmd + '=' + str(toks.split(',')))

        try:
            exec('glbl.' + cmd + ' = ' + arg)
        except:
            exec('glbl.' + cmd + ' = ' + str(arg.split(',')))

    #==============================================================================
    # add Parameter to the Parameters class instance
    #==============================================================================
    def addGlobalParameters(self, runNumber):

        # global errors,  global_parms, parms

        for parmName in self.global_parms:
            oldparm = self.parms[parmName + '_1']
            self.parms.add(parmName + '_' + str(runNumber),
                      value = oldparm.value,
                      vary = oldparm.vary,
                      min = oldparm.min,
                      max = oldparm.max,
                      expr = parmName + '_1')


    #==============================================================================
    # process end of section
    #==============================================================================
    def processEndSection(self, runNumber, lineNumber, glbl):

        optionalList = ['n_delt']
        
        # check if a signal file was specified
        if not glbl.signal_filename:
            print('\n *** Error: No signal filename for')
            print(' *** data section ending at line ', lineNumber)
            self.errors.append('no signal filename for section ending at line ' 
                + str(lineNumber))
            

        # if this is not the first section, add global parameters
        #   and duplicate old list items if needed
        if runNumber > 1:
            self.addGlobalParameters(runNumber)

            # if var or tuple item was not supplied, use value from previous run if available
            for item in self.lists_list:
                if not eval('glbl.' + item):
                    try:
                        exec('glbl.' + item + ' = glbl.' + item + 's[' + str(runNumber -2) + ']')
                    except:
                        pass
            for item in self.tuples_list:
                if not eval('glbl.' + item):
                    try:
                        exec('glbl.' + item + ' = glbl.' + item + 's[' + str(runNumber -2) + ']')
                    except:
                        pass

        # add item to items list and check if required items have been supplied
        for item in self.lists_list:
            exec('glbl.' + item + 's.append(glbl.' + item + ')')
            if not eval('glbl.' + item) and not item in optionalList:
                print('\n *** Error: No ', item, ' for')
                print(' *** data section ending at line ', lineNumber)
                self.errors.append('no ' + str(item) + ' for section ending at line ' +str(lineNumber))

        # add tuple to tuples list
        for item in self.tuples_list:
            exec('glbl.' + item + 's.append(glbl.' + item + ')')

        # make a list of required parameters and then check if they have been supplied
        required_parms = ['y_scale', 'baseline', 'ion_tof', 'ffr', 'temp']

        if glbl.cutoff_type:        
            if glbl.cutoff_type == 'exp':
                required_parms.append('ecutm')
                required_parms.append('ecuts')
    
            elif glbl.cutoff_type == 'tanh':
                required_parms.append('tcutc')
                required_parms.append('tcutw')
        
        if glbl.function:
            if glbl.function == 'erf':
                required_parms.append('e0')
                required_parms.append('w')


        # test if required parameters have been supplied
        # if parameter is missing, use parameter from previous run number
        for p in required_parms:
            try:
                self.parms[p + '_' + str(runNumber) ].value
                continue
            except:
                if runNumber > 1:
                    oldparm = self.parms[p + '_' + str(runNumber - 1)]
                    self.parms.add(p + '_' + str(runNumber),
                      value = oldparm.value,
                      vary  = oldparm.vary,
                      min   = oldparm.min,
                      max   = oldparm.max,
                      expr  = oldparm.expr)

                else:
                    print('\n *** Error: No parameter', p, ' for')
                    print(' *** data section ending at line ', lineNumber)
                    self.errors.append('no parameter ' + p + ' for section ending at line ' +str(lineNumber))

        # save filenames and reset to None for entry to next section
        glbl.signal_filenames.append(glbl.signal_filename)
        glbl.background_filenames.append(glbl.background_filename)
        glbl.signal = None
        glbl.background = None
        glbl.function = None
        glbl.averaging_type = None
        glbl.cutoff_type = None
        glbl.Tmin = None
        glbl.Tmax = None




    #==============================================================================
    # parse command file
    #==============================================================================
    def parse_cmd_file(self, filename, glbl, parms):
        self.parms = parms
        runNumber = 1

        # open and read command file
        cmdFile=self.openForReading(filename)
        lines = cmdFile.readlines()
        cmdFile.close()

        # parse command file, line by line
        for i in range(len(lines)):
            line=lines[i].strip()
            lineNumber = i+1

            # skip blank lines and comments
            if line.startswith('#') or len(line) == 0:
                continue

            # break line into command, args, and tokens of the args
            cmd, arg, tokens = self.getTokens(line)

            # check if line specifies a parameter
            if self.isParameter(cmd):
                self.addParameter(runNumber, cmd, arg, tokens, line, lineNumber)

            # check if line specifies a variable
            elif self.isVar(cmd):
                self.addVar(cmd, arg, tokens, line, lineNumber, glbl)

            # check if line specifies an item to be added to a list of values
            #   that vary with run number
            elif self.isListItem(cmd):
                self.addListItem(cmd, arg, tokens, line, lineNumber, glbl)

            # check if line specifies a tuple to be added to a list (e.g. (t_min, t_max))
            elif self.isTupleItem(cmd):
                self.addTupleItem(cmd, arg, tokens, line, lineNumber, glbl)

            # check if line specifies a signal or background file name
            # note this could be consolidatae with lists_list items
            elif cmd == 'signal' or cmd == 'background':
                #use following code to input file name without quotes
                filelist = []
                for filename in glob(tokens[0]):
                    filelist.append(filename)

                if len(filelist) == 0:
                    print('***** error no file found \n')
                    print('*****', tokens[0])
                    self.errors.append(' no match for file patern' + tokens[0] +
                                    ' on line ' + str(lineNumber))
                    filename = 'no match for pattern'

                if len(filelist) > 1:
                    print('***** error nonunique match to file pattern \n')
                    print('*****', tokens[0])
                    print('***** filelist =', filelist)
                    self.errors.append(' nonunique match for file pattern' + tokens[0] +
                                    ' on line ' + str(lineNumber))
                    return
                exec ('glbl.' + cmd + '_filename =' + "filename")

            # check if line indicates end of dataset section
            elif cmd.upper().startswith('END'):
                self.processEndSection(runNumber, lineNumber, glbl)
                runNumber += 1

            else:
                print('\n *** unknown command on line', i+1)
                print(' *** command =', cmd)
                self.errors.append('unknown command on line ' + str(i+1) + ' command = ' + cmd)

        return self.errors


#==============================================================================
# test parseCmdFile if executed as main
#==============================================================================
if __name__ == '__main__':         
    
    
    # global errors,  global_parms, parms

    glbl = TOF_fit_global.TOF_fit_global()
    parms = Parameters2()
    fit_control = Fit_control()

    print()
    print(' Test ParseCmdFile')
    print()
    print('=====================')
    print('    Initial State    ')
    print('=====================')
    fit_control.print_list_items()
    fit_control.print_tuple_items()

#==============================================================================
#     # parse the constants file    
#     const_filename = 'testSetVariables.dat'   
#     parms, functions, signal_filename, backgrounds, errors = parseCmdFile(const_filename)
#==============================================================================

    filename = 'Fits\\fit0011.fit_in'

    errors = fit_control.parse_cmd_file(filename, glbl, parms)
    
    print()
    print('=========================')    
    print('       Final State       ')
    print('=========================')
    print('signal_filenames:')
    for file in glbl.signal_filenames:
        print('  ', file)

    print()    
    print('background_filenames:')
    for file in glbl.background_filenames:
        print('  ', file)

    print()
    print('fit_ranges:      ', glbl.fit_ranges)
    print('comment_xlsx     ', glbl.comment_xlsx)
    print()
    
    fit_control.print_list_items()
    fit_control.print_tuple_items()
    print()
    
    is_global_fit = any([parms[p].glbl for p in parms])
    print('{:16s} = {}'.format('Global Fit', is_global_fit))
    print()
    
    parms.pretty_print()
    
    if len(errors) > 0:
        print('*****')
        print('***** Errors in command file found')
        for error in errors:
            print('***** ', error)
        print('*****')
        
        raise SystemExit ('Error in command file')
    
    