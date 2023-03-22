delete all 

load /media/rhys/data11/ampk_rhys_2022/metadynamics/data/a2b1+A769/01-Min/a2b1+A769_minimised.pdb, a2b1
load /media/rhys/data11/ampk_rhys_2022/metadynamics/data/a2b2+A769/01-Min/a2b2+A769_minimised.pdb, a2b2

align a2b1 and resi 5-360 and name ca, a2b2 and resi 5-360 and name ca

select s0_1, a2b1 and (resi 81 and name ca or resi 11 and name ca) 
select s1_1, a2b1 and (resi 302 and name ca or resi 43 and name ca) 
select s0_2, a2b2 and (resi 81 and name ca or resi 11 and name ca) 
select s1_2, a2b2 and (resi 303 and name ca or resi 43 and name ca)

run ~/pymol_scripts/com.py

com s0_1, object=p0_1
com s1_1, object=p1_1
com s0_2, object=p0_2
com s1_2, object=p1_2

run group_tools/PyMOL_scripts/draw_funnel.py

draw_funnel p0_1, p1=p0_1, p2=p1_1, s_cent=25, beta_cent=0.10, wall_width=15, wall_buffer=1.5, lower_wall=0, upper_wall=45, vec_step=2.5, angle_sample=18
set_name funnel, funnel_a2b1
draw_funnel p0_2, p1=p0_2, p2=p1_2, s_cent=25, beta_cent=0.10, wall_width=15, wall_buffer=1.5, lower_wall=0, upper_wall=45, vec_step=2.5, angle_sample=18
set_name funnel, funnel_a2b2


#hide lines, all
#show cartoon, pdb1
#show sticks, pdb1 and resn FRG
#show cartoon, pdb2
#show sticks, pdb2 and resn FRG
#set sphere_scale, 0.4
#show spheres, funnel
#
#color white, pdb1
#color lightorange, pdb2
#color magenta, (resn FRG)
#color forest, funnel


