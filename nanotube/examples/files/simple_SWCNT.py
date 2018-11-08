from nanotube import SWCNT

armchair = SWCNT(3,radius=.5)
armchair.save('armchair.mol2')

zigzag = SWCNT(3,radius=.5,chirality='zigzag')
zigzag.save('zigzag.mol2')

armchair2 = SWCNT(4,n=8,m=8)
armchair2.save('armchair2.mol2')