import mbuild as mb
import numpy as np
from copy import deepcopy
from math import *

__all__ = ['SWCNT']

class SWCNT(mb.Compound):
""" A single-walled Carbon nantube recipe

n,m : Chriality parameters
"""
    def __init__(self,n,m):
        super(SWCNT, self).__init__()

        class C(mb.Compound):
            def __init__(self):
                super(C,self).__init__()
                
                self.add(mb.Particle(name='C'))

        cell = mb.Compound()

        unit_length = .071

        if n != m:
            for i in range(0,n*2):
                c = C()
                cell.add(c)
                if i % 2 == 0:
                    c.translate((sqrt(3)*unit_length*i,0,0))
                else:
                    c.translate((sqrt(3)*unit_length*i,unit_length,0))