import mbuild as mb
import numpy as np
from copy import deepcopy
from math import *

__all__ = ['SWCNT']

class SWCNT(mb.Compound):
    """ A single-walled Carbon nanotube recipe
    radius: radius of nanotube in nm
    chirality : If using radius instead of n and m, the type of chirality desired
        either "armchair" or "zigzag". "armchair" by defualt

    n,m : chirality parameters, if only n is defined armchair is assumed
    length : length of nanotube in nm

    note: either the radius or the chiraity parameters must be defined but not
    both.
    """

    def __init__(self,length,radius=None,chirality="armchair",n=None,m=None):
        super(SWCNT, self).__init__()

        nanotube = mb.Compound()

        # define repeatable row units for the different chiralities
        _ARMCHAIR = np.array([(0,0,0),
                    (.071,.071*np.sqrt(3),0),
                    (.213,.071*np.sqrt(3),0),
                    (.284,0,0)])

        _ARM_WALK = .426

        _ZIGZAG = np.array([(0,0,0),
                    (.071*np.sqrt(3),.071,0)])

        _ZIG_WALK = 2*.071*np.sqrt(3)

        _CHIRAL_UNIT = {'armchair':_ARMCHAIR, 'zigzag':_ZIGZAG}
        _CHIRAL_LENGTH = {'armchair':_ARM_WALK, 'zigzag':_ZIG_WALK}
        _CHIRAL_HEIGHT = {'armchair': .246, 'zigzag': .142}

        _row = []
        _sheet = []

        # define Carbon atom
        class C(mb.Compound):
            def __init__(self,position):
                super(C,self).__init__(pos=position)
                
                self.add(mb.Particle(name='C'))

        #Checking validity of input values
        if radius is None and n is None:
            raise ValueError("Either the radius or the chiraity parameters"
            "must be defined")
        elif radius is not None and n is not None:
            raise RuntimeWarning("Both radius and chirality parameters defined,"
            "radius overriding (n,m) values")

            n = None
            m = None

        elif radius is not None and radius < 0.19:
            raise ValueError("The smallest possible radius is 0.19 nm. Your"
            " radius was " + str(radius) + " nm.")
        elif m is not None and not(m is not n or m is not 0):
            raise ValueError("m must equal n or 0. Chiral nanotubes not "
            "currently supported")

        if radius is not None and not(chirality in ['armchair','zigzag']):
            raise ValueError("Chirality must be either \'armchair\' or"
            "\'zigzag\'")


        #Sets implied chirality from the n,m values if nessecary
        if n is not None and m is 0:
            chirality = "zigzag"

        #Find number of cells in row and number of rows in tube
        if radius is None:
            n_cells = n
        else:
            #Divides circumference of the desired nanotube by the space taken
            # up by one repeatable unit for specified chirality  
            n_cells_aprx = nano_radius*(np.pi*2)/_CHIRAL_LENGTH[chirality]

            n_cells = int(round(n_cells_aprx))

        #Finds number of rows in tube
        z = int(round((length-.071)/.213) + 1) # WIP

        #Propagate out to a row
        for i in range(0,n):
            temp = deepcopy(_CHIRAL_UNIT[chirality])
            for Carbon in temp:
                Carbon[0] = Carbon[0]+_CHIRAL_LENGTH[chirality]*i
                _row.append(Carbon)

        #Propagate out to a sheet

        for i in range(0,z): 
            temp =  deepcopy(_row)
            for Carbon in temp:
                if chirality == 'zigzag' and i%2 == 0 and (Carbon[0] % _ZIG_WALK) < 0.01:
                    Carbon[1] = Carbon[1] *-1
                Carbon[1] = Carbon[1] + _CHIRAL_HEIGHT[chirality]*i
                
            _sheet.append(temp)

        for row in _sheet:
            for Carbon in row:
                atom = C(Carbon)
                nanotube.add(atom)

        nanotube.save("nanotube.mol2")

        # -> folded tubes 
        # -> insert C atoms
        # 
        # Afterwards take care of eccentricity factor for small tubes (< 14 C atoms per circle)
        