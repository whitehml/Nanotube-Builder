from nanotube import SWCNT_solvated
import mbuild as mb

water = mb.load('files/tip3p.mol2')
water.name = 'water'

slow = SWCNT_solvated(solv=water,density=.42,radius=1.5)
slow.save('solvated_SWCNT_slow.mol2')

quick = SWCNT_solvated(solv=water,density=.42,radius=1.5,superPacked=True)
quick.save('solvated_SWCNT_quick.mol2')