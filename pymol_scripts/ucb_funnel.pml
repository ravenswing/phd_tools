#load /workdrive/Project/ucb_development_project/BindingTrajectory/5ai5/input.pdb, pdb1
delete all 

load /media/rhys/ExtHD/Project/ucb_development_project/BindingTrajectory/first_pocket.pdb, pdb1
load /media/rhys/ExtHD/Project/ucb_development_project/BindingTrajectory/second_pocket.pdb, pdb2

align pdb1 and resi 226-546, pdb2 and resi 226-546 

# centred between pockets
#select p2, resi 335 and name cg
#select p1, resi 243 and name ca

select p1, resi 238 and name ca
select p2, resi 270 and name ca

run group_tools/PyMOL_scripts/draw_funnel.py

#draw_funnel p2, p1, p2, wall_width=22, upper_wall=40, s_cent=26
draw_funnel p2, p1, p2, wall_width=23, lower_wall=5, upper_wall=55, s_cent=38

hide lines, all
show cartoon, pdb1
show sticks, pdb1 and resn FRG
show cartoon, pdb2
show sticks, pdb2 and resn FRG
set sphere_scale, 0.4
show spheres, funnel

color white, pdb1
color lightorange, pdb2
color magenta, (resn FRG)
color forest, funnel


