import mbuild as mb
import numpy as np
from copy import deepcopy
from math import *

__all__ = ['SWCNT']

class SWCNT(mb.Compound):
""" A single-walled Carbon nantube recipe

n,m : Chriality parameters
length : length of nanotube in nm
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
                    
        cell_2 = deepcopy(cell)
        mb.rotate(cell_2, pi, around=np.asarray([1, 0, 0]))
        cell_2.translate((0,-.142,0))

        unit_cell = mb.Compound()
        unit_cell.add(cell)
        unit_cell.add(cell_2)
        unit_cell.translate((0,.213,0))

        #Fold the unit cell
        circumference = .246 * n

        for C in unit_cell:
            x_index = (C.pos[0]/circumference)*2*pi
            x = cos(x_index) * (circumference/(2*pi))
            y = sin(x_index) * (circumference/(2*pi))
            z = C.pos[1]
            
            C.translate_to((x,y,z))

        tube = mb.Compound()
        tube.add(unit_cell)

        chain_length = int(length/.284)
        for i in range(1,chain_length):
            prev = tube.children[-1]
            new_cell = mb.clone(unit_cell)
            if i % 2 == 1:
                new_cell.rotate(pi, around=np.asarray([0, 1, 0]))
                new_cell.translate_to((0,0,prev.pos[2] + .213))
            else:
                new_cell.translate((0,0,prev.pos[2] - .03))
            tube.add(new_cell)