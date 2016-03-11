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
        self.parmList = ['A', 'E0', 'W',
                    'TCutC', 'TCutW', 'ECutM', 'ECutS',
                    'Temp', 'FFR', 'IonTOF', 'Yscale', 'Baseline']

        self.constList =['ThetaStep', 'ZDiffWall', 'RDiffWall',
                    'ZRef', 'ZAperture', 'RAperture',
                    'ZSource', 'RSource', 'ZLaser', 'LLaser', 'ZFinal', 'RFinal',
                    'NPointsDetector', 'NPointsSource', 'GridType',
                    'comment_xlsx']

        self.listList = ['DataColLine', 'DataRowLine', 'MoleculeLine', 'TemperatureLine',
                    'VibStateLine', 'RotStateLine',
                    'Tmin', 'Tmax', 'Function', 'AveragingType', 'cutoff_type']

        self.errors = []
        self.globalParms=[]
        self.parms = None

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

        except:
            quit


    #==============================================================================
    # break input line into tokens
    #=========================================================================
    def getTokens(self, line):

        # split the line on ','
        tokens = line.split(',')
        for j in range(len(tokens)):
        # replace nule tokens (i.e. '') with None
            tokens[j] = tokens[j].strip()
            if tokens[j] == '': tokens[j] = None

            try:
                #test if token has an embedded '#'
                ii = tokens[j].index('#')
                # if so, remove rest of token and strip white space
                tokens[j] = tokens[j][0:ii].strip()
                # ignore rest of tokens
                for k in range(j + 1, len(tokens)):
                    tokens.pop(-1)
                break
            except:
                pass
        return tokens


    #==============================================================================
    # check if cmd is an add Parameter item
    #==============================================================================
    def isParameter(self, tokens):
        return any(tokens[0] in s for s in self.parmList)


    #==============================================================================
    # check if cmd is an add Global Variable item
    #==============================================================================
    def isGlobalConst(self, tokens):
        matching = [tokens[0] == s for s in self.constList]
        return any(matching)


    #==============================================================================
    # check if cmd is a List item
    #==============================================================================
    def isListItem(self, tokens):
        matching = [tokens[0].upper() == s.upper() for s in self.listList]
        #print('isListItem: tokens=',tokens, 'any=', any(matching)) #, 'matching=', matching)
        return any(matching)

    #==============================================================================
    # check if token is a number
    #==============================================================================
    def is_number(self, s):
        try:
            float(s)
            return True
        except:
            return False

    #==============================================================================
    # add Parameter to the Parameters class instance
    #==============================================================================
    def addParameter(self, runNumber, tokens, line, lineNumber):

        # global errors,  globalParms, parms

        # check if number of tokesn is correct
        #   there should be name, and 2, 3 or 4 parameters
        #   i.e. 3 to 5 tokens on the line
        if len(tokens) < 3 or len(tokens) > 5:
            print('\n *** wrong number of items on line ', lineNumber)
            print(' *** items found =', len(tokens),
                  '. Should be 3, 4, or 5')
            print(' *** line =', line)
            print(' *** items =', tokens)
            self.errors.append(('wrong number of items', 'line ' + str(lineNumber)))
            return

        # set defaults
        vary1 = False
        glbl1 = False
        min1 = None
        max1 = None

        # set vary and global
        if tokens[2] != '0':
                vary1 = True
        if tokens[2] == '2':
            glbl1 = True
            self.globalParms.append(tokens[0])

        # set bounds
        try:
            min1=float(tokens[3])
        except:
            pass
        try:
            max1 = float(tokens[4])
        except:
            pass

        self.parms.add(tokens[0] + '_' + str(runNumber),
                  value = float(tokens[1]),
                  vary = vary1,
                  glbl = glbl1,
                  min = min1,
                  max = max1
                  )
        # print('addParameter runNumber, tokens =', runNumber, tokens)


    #==============================================================================
    # add const
    #==============================================================================
    def addGlobalConst(self, tokens, line, lineNumber, glbl):
        if len(tokens) != 2:
            error = 'wrong number of items for const on line ' + str(lineNumber)
            print('\n ***' + error, lineNumber)
            print(' *** items found =', len(tokens),
                  '. Should be 2')
            print(' *** line =', line)
            print(' *** items =', tokens)
            self.errors.append(error)
            return

        exec('glbl.' + tokens[0] + '=' + tokens[1])


    #==============================================================================
    # add list item
    #==============================================================================
    def addListItem(self, tokens, line, lineNumber, glbl):

        # global errors,  globalParms, parms

        if len(tokens) != 2:
            error = 'wrong number of arguments items on line ' + str(lineNumber)
            self.errors.append(error)

            print('\n ***' + error)
            print(' *** items found =', len(tokens),
                  '. Should be 2')
            print(' *** line =', line)
            print(' *** items =', tokens)


        # item is a string, need to add quotes to exec, e.g. ERF become 'ERF'
        if self.is_number(tokens[1]):
            exec('glbl.' + tokens[0] + '=' + tokens[1])
        else:
            exec('glbl.' + tokens[0] + '=' + "tokens[1]")

    #==============================================================================
    #     # append the item to list
    #         cmd = 'glbl.' + tokens[0] + 's.append(glbl.' + tokens[0] +')'
    #         exec(cmd)
    #         # exec('glbl.' + tokens[0] + 's.append(glbl.' + tokens[0])
    #         pass
    #==============================================================================

    #==============================================================================
    # add Parameter to the Parameters class instance
    #==============================================================================
    def addGlobalParameters(self, runNumber):

        # global errors,  globalParms, parms

        for parmName in self.globalParms:
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
        # global errors,  globalParms, parms

        optionalList = ['Tmin', 'Tmax']

        # if this is not the first section, add global parameters
        #   and duplicate old list items if needed
        if runNumber > 1:
            self.addGlobalParameters(runNumber)

            # if item was not supplied, use value from previous run if available
            for item in self.listList:
                if not eval('glbl.' + item):
                    try:
                        exec('glbl.' + item + ' = glbl.' + item + 's[' + str(runNumber -2) + ']')
                    except:
                        pass




        # add item to items list and check if required items have been supplied
        for item in self.listList:
            exec('glbl.' + item + 's.append(glbl.' + item + ')')
            if not eval('glbl.' + item) and not item in optionalList:
                print('\n *** Error: No ', item, ' for')
                print(' *** data section ending at line ', lineNumber)
                self.errors.append('no ' + str(item) + ' for section ending at line ' +str(lineNumber))

        # check if required parameters have been supplied
        required_parms = ['Yscale', 'Baseline', 'IonTOF', 'FFR', 'Temp']

        if glbl.cutoff_type.lower() == 'exp':
            required_parms.append('ECutM')
            required_parms.append('ECutS')

        elif glbl.cutoff_type.lower() == 'tanh':
            required_parms.append('TCutC')
            required_parms.append('TCutW')

        if glbl.Function.lower() == 'erf':
            required_parms.append('E0')
            required_parms.append('W')


        for p in required_parms:
            try:
                self.parms[p + '_' + str(runNumber) ].value
                continue
            except:
                if runNumber > 1:
                    oldparm = self.parms[p + '_' + str(runNumber - 1)]
                    self.parms.add(p + '_' + str(runNumber),
                      value = oldparm.value,
                      vary = oldparm.vary,
                      min = oldparm.min,
                      max = oldparm.max,
                      expr = oldparm.expr)
                else:
                    print('\n *** Error: No parameter', p, ' for')
                    print(' *** data section ending at line ', lineNumber)
                    self.errors.append('no parameter' + p + ' for section ending at line ' +str(lineNumber))


        # save filenames and reset to None for entry to next section
        glbl.signalFiles.append(glbl.signalFile)
        glbl.backgroundFiles.append(glbl.backgroundFile)
        glbl.signalFile = None
        glbl.backgroundFile = None
        glbl.Function = None
        glbl.AveragingType = None
        glbl.cutoff_type = None
        glbl.Tmin = None
        glbl.Tmax = None




    #==============================================================================
    # parse command file
    #==============================================================================
    def parse_cmd_file(self, filename, glbl, parms):

        # global errors,  globalParms, parms

        self.parms = parms
        # initialize()
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

            # break line into tokens
            tokens = self.getTokens(line)

            # check if line specifies a parameter
            if self.isParameter(tokens):
                self.addParameter(runNumber, tokens, line, lineNumber)

            # check if line specifies a constant
            elif self.isGlobalConst(tokens):
                self.addGlobalConst(tokens, line, lineNumber, glbl)

            elif self.isListItem(tokens):
                self.addListItem(tokens, line, lineNumber, glbl)
                pass

            # check if line specifies a signal file name
            # note this could be consolidatae with listList items
            elif tokens[0].upper() == 'SIGNAL':
                #use following code to input signal file name without quotes
                filelist = []
                for file in glob(tokens[1]):
                    filelist.append(file)

                if len(filelist) == 0:
                    print('***** error no file found \n')
                    print('*****', tokens[1])
                    self.errors.append(' no match for file patern' + tokens[1] +
                                    ' on line ' + str(lineNumber))
                    file = 'no match for pattern'

                if len(filelist) > 1:
                    print('***** error nonunique match to file pattern \n')
                    print('*****', tokens[1])
                    print('***** filelist =', filelist)
                    self.errors.append(' nonunique match for file patern' + tokens[1] +
                                    ' on line ' + str(lineNumber))
                    return

                glbl.signalFile = file

            # check if line specifies a background file name
            # note this could be consolidatae with listList items
            elif tokens[0].upper() == 'BACKGROUND':
                # use following code to input background file name without quotes
                filelist = []
                for file in glob(tokens[1]):
                    filelist.append(file)
                if len(filelist) > 1:
                    print('***** error nonunique match to file pattern \n')
                    print('*****', tokens[1])
                    print('***** filelist =', filelist)
                    self.errors.append(' nonunique match for file patern' + tokens[1] +
                                    ' on line ' + str(lineNumber))

                glbl.backgroundFile = file
                # use the following line to file name with quotes
                # exec('glbl.backgroundFile =' + tokens[1])
                continue

            # check if line indicates end of dataset section
            elif tokens[0].upper().startswith('END'):
                self.processEndSection(runNumber, lineNumber, glbl)
                runNumber += 1

            else:
                print('\n *** unknown command on line', i+1)
                print(' *** command =', tokens[0])
                self.errors.append('unknown command on line ' + str(i+1) + ' command = ' + tokens[0])

        return self.errors


    #--------------------------------------------------------------------------------------------------
    # print list items
    #--------------------------------------------------------------------------------------------------
    def print_list_items(self):
        for item in self.listList:
            # print('{:16s} = {}'.format(item + 's', eval('glbl_old.' + item + 's')))
            print('{:16s} = {}'.format(item + 's', eval('glbl.' + item + 's')))



    
#==============================================================================
# test parseCmdFile if executed as main
#==============================================================================
if __name__ == '__main__':         
    
    
    # global errors,  globalParms, parms

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
    print()    

#==============================================================================
#     # parse the constants file    
#     const_filename = 'testSetVariables.dat'   
#     parms, functions, signalFiles, backgroundFiles, errors = parseCmdFile(const_filename)
#==============================================================================

    filename = 'Fits\\fit031.tof_in'



    # parse_result = parseCmdFile(filename)
    
    errors = fit_control.parse_cmd_file(filename, glbl, parms)
    
    print()
    print('=========================')    
    print('       Final State       ')
    print('=========================')
    print()
    print('glbl.Functions:           ', glbl.Functions)
    
    print('glbl.signalFiles:         ')
    for file in glbl.signalFiles:
        print('  ', file)

    print('glbl.backgroundFiles:     ')
    for file in glbl.backgroundFiles:
        print('  ', file)

    print('glbl.Tmins:               ', glbl.Tmins)
    print('glbl.Tmaxs:               ', glbl.Tmaxs)
    print('glbl.comment_xlsx         ', glbl.comment_xlsx)
    print()
    
    
    fit_control.print_list_items()
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
    
    