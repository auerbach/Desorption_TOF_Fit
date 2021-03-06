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
                          'energy_angle_scaling',
                          'grid_type',
                          'length_laser',
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
                          'z_ref',
                          'z_source'
                          ]

        self.lists_list = ['averaging_type',
                           'comment',
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
    # print variable items
    #==========================================================================
    def print_var_items(self):
        for item in self.vars_list:
            print('{:16s} = {}'.format(item + 's', eval('glbl.' + item)))
    
    #==========================================================================
    # print list items
    #==========================================================================
    def print_list_items(self):
        for item in self.lists_list:
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


    #==========================================================================
    # break string on seperator and remove imbedded comments
    #=========================================================================
    def tokenize(self, line, sep, count):
        toks = line.split(sep, count)
        for j in range(len(toks)):
            toks[j] = toks[j].strip()
            
            # replace null tokens (i.e. '') with None
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

        
    #==========================================================================
    # use tokenize to get tokens from input line
    #
    # first split line on '='; then split remainder on ','
    #=========================================================================    
    def getTokens(self, line):
        # first split the line on '=,'
        toks = self.tokenize(line, '=', 1)
        cmd = toks[0].lower()
        
        # if there is more than 1 token break apart further
        # (sometimes there will only be one as in 'End Section')
        if len(toks)>1:
            arg = toks[1].lower()
            # split the args on ',' to get tokens
            tokens = self.tokenize(arg, ',', -1)
        else:
            arg = None
            tokens = None

        return cmd, arg, tokens


    #==============================================================================
    # check if tokens specify a Parameter item
    #==============================================================================
    # def isParameter(self, tokens):
    #     return any(tokens[0] in s for s in self.parms_list)
    def isParameter(self, cmd):
        return any(cmd in s for s in self.parms_list)


    #==============================================================================
    # check if tokens specify a Global Variable item
    #==============================================================================
    def isVar(self, cmd):
        return any([cmd == s for s in self.vars_list])


    #==============================================================================
    # check if tokens specify a List item
    #==============================================================================
    def isListItem(self, cmd):
        matching = [cmd == s for s in self.lists_list]
        return any(matching)

    #==============================================================================
    # check tokens specify a Tuple item
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
        # check for correct number of parameters
        #   cmd is the parameter name; tokens are the parameters
        #   there should be name, and 2, 3 or 4 parameters
        if len(tokens) < 2 or len(tokens) > 4:
            print('\n *** wrong number of values for Parameter Item on line ', lineNumber)
            print(' *** number found =', len(tokens),
                  '. Should be 2, 3, or 4')
            print(' *** line =', line)
            print(' *** items =', tokens)
            self.errors.append(('wrong number of items', 'line ' + str(lineNumber)))
            return

        # set defaults whether parameter varies or is fixed; is global or not;
        #   and whether there is a min and max constraint
        vary1 = False
        glbl1 = False
        min1 = None
        max1 = None

        # token 2 to set vary and global;
        #   0 -> fixed; 1 -> vary; 2-> vary and global
        if tokens[1] != '0':
                vary1 = True
        if tokens[1] == '2':
            glbl1 = True
            # store the name of global parameters in a list.  This could,
            #   and probably should be replaced by using the glbl property 
            #   of parameters to process global parameters
            #   keep for now because it is simple and works
            self.global_parms.append(cmd)

        # set bounds if min or max is given
        try:
            min1=float(tokens[2])
        except:
            pass
        try:
            max1 = float(tokens[3])
        except:
            pass

        # add the parameter to the parameter dictionary        
        self.parms.add(cmd + '_' + str(runNumber),
                  value = float(tokens[0]),
                  vary = vary1,
                  glbl = glbl1,
                  min = min1,
                  max = max1
                  )


    #==============================================================================
    # add variable
    #==============================================================================
    def addVar(self, cmd, arg, tokens, line, lineNumber, glbl):
        # check there is one value to use to set the variable
        if len(tokens) != 1:
            error = 'wrong number of args for item on line ' + str(lineNumber)
            print('\n ***' + error, lineNumber)
            print(' *** items found =', len(tokens),
                  '. Should be 1')
            print(' *** line =', line)
            print(' *** items =', tokens)
            self.errors.append(error)
            return
            
        # use python exec command to execute line var = value
        # exec('glbl.' + cmd + '=' + tokens[0])
        if self.isNumber(tokens[0]):
            exec('glbl.' + cmd + '=' + tokens[0])
        else:
            exec('glbl.' + cmd + '=' + "tokens[0]")


    #==============================================================================
    # add list item
    #==============================================================================
    def addListItem(self, cmd, arg, tokens, line, lineNumber, glbl):

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
    # duplicate global parameters for run >
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

        optionalList = ['n_delt', 'comment']
        no_dup_list  = ['comment']
        
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
                if item in no_dup_list:
                    continue
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
        required_parms = ['y_scale', 'baseline', 'ion_tof', 'ffr']

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
        glbl.comment = None




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
                    self.errors.append(' no match for file patern ' + tokens[0] +
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
    
    import subprocess
    
    # -------------------------------------------------------------------------
    # create instances of TOF-fit_global, Parameters2 and Fit_control classes
    # -------------------------------------------------------------------------
    glbl = TOF_fit_global.TOF_fit_global()
    parms = Parameters2()
    fit_control = Fit_control()

    print()
    print(' Test ParseCmdFile')
    print()
    print('=====================')
    print('    Initial State    ')
    print('=====================')
    fit_control.print_var_items()
    fit_control.print_list_items()
    fit_control.print_tuple_items()

    #==============================================================================
    #   get path to the command file
    #==============================================================================
    with open("path_to_fits_and_editor.txt") as file: 
        lines = file.readlines()
        path_to_fits = lines[0].split(':',1)[1].strip()
        editor_cmd   = lines[1].split(':',1)[1].strip()
        
    fit_number = '0036'
    prompt = True  
    #------------------------------------------------------------------------------
    # begin1 comment out for testing
    #------------------------------------------------------------------------------
    
    if prompt:
        # Get Last fit number from FitNumber.dat
        fn = path_to_fits + 'FitNumber.dat'
        print('fn =', fn)
        with open(fn, 'r+') as fit_number_file:
            fit_number = '{:04d}'.format(int(fit_number_file.readline()))
            prompt = 'please enter fit number to test parsing: [' + fit_number + ']'
            
            while True:
                ans = input(prompt)
                if ans:
                    n = '{:04d}'.format(int(ans))
                else:
                    n = fit_number
            
                if int(n) > int(fit_number):
                    print('maximum available fit number is ', fit_number)
                else:
                    break
            
        fit_number = n
        cmd_fn = path_to_fits + 'Fit' + fit_number + '.fit_in'
        
        ans = input('edit file? [no]')
        if ans.lower().startswith('y'):   
            subprocess.run(editor_cmd + ' "' + cmd_fn + '"')

    errors = fit_control.parse_cmd_file(cmd_fn, glbl, parms)
    
    print()
    print('=========================')    
    print('       Final State       ')
    print('=========================')
    fit_control.print_var_items()
    fit_control.print_list_items()
    fit_control.print_tuple_items()
    print()

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
    print('comments         ', glbl.comments)
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
    
    