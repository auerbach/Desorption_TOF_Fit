# -*- coding: utf-8 -*-

from lmfit import Parameter, Parameters


###############################################################################
# class Parameters -- derived from lmfit.Parameter
###############################################################################
class Parameter2 (Parameter):
    """
    derived from lmfit.Parmeter
        adds boolen attribute glbl_old to specify if parameter is global or not
        adds support for glbl_old in pretty_print()
    """
    
    def __init__(self, glbl = False, *args, **kwargs):
        Parameter.__init__(self, *args, **kwargs)
        #super().__init__(self, *args, **kwargs)
        self.glbl = glbl
        
    def __repr__(self):
        s = []
        if self.name is not None:
            s.append("'%s'" % self.name)
            
        sval = repr(self._getval())
        
        if self.stderr is not None:
            sval = "value=%s +/- %.3g" % (sval, self.stderr)
        
        elif not self.vary and self._expr is None:
            sval = "value=%s (fixed" % (sval)
            sbound = None
        else:
            sval = "value=%s (varys" % (sval)
            sbound = "bounds=[%s:%s]" % (repr(self.min), repr(self.max))
        s.append(sval)
        
        #if(self.glbl_old):
        #sglbl = "self.glbl_old=%s" % (self.glbl_old)
        s.append("global =%s" %(self.glbl))
        
        if sbound:
            s.append(sbound)
        
        
        if self._expr is not None:
            s.append("expr='%s'" % (self.expr))
        return "<Parameter %s>" % ', '.join(s)



###############################################################################
# class Parameters2 -- derived from lmfit.Parameters
###############################################################################
class Parameters2(Parameters):
    """
    Parameters2 is dervied from lmfit.Parameters
        overrides the add() method to support global attribute of Parameter2
        
    """
        
    #def __init__(self, asteval=None, *args, **kwargs):
    #    Parameters.__init__(asteval, *args, **kwargs)
        
    def add(self, name, value=None, vary=True, glbl=False, min=None, max=None, expr=None):
        """
        Convenience function for adding a Parameter:

        Example
        -------
        p = Parameters()
        p.add(name, value=XX, glbl_old=True ...)

        is equivalent to:
        p[name] = Parameter(name=name, value=XX, glbl_old=True....
        """
        if isinstance(name, Parameter2):
            self.__setitem__(name.name, name)
        else:
            self.__setitem__(name, Parameter2(value=value, name=name, vary=vary, glbl=glbl,
                                             min=min, max=max, expr=expr))


if __name__ == '__main__':

    # test original Parameters class
    parms = Parameters()
    parms.add('p1', value = 1, vary = True)
    parms.add('p2', value = 2, vary = False)
    parms.pretty_print()
    
    # test new Parameters class with global attriute
    parms2 = Parameters2()
    parms2.add('p1', value = 1, vary = True, glbl=True)
    parms2.add('p2', value = 2, vary = False, glbl=False)
    
    for p in parms2:
        print(parms2[p].name, parms2[p].vary, parms2[p].glbl)
        
    parms2.pretty_print()



