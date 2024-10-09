#delete all 
#run ~/pymol_scripts/com.py
#run group_tools/PyMOL_scripts/draw_funnel.py

#load /media/rhys/data11/ampk_rhys_2022/metadynamics/data/a2b1+A769/01-Min/a2b1+A769_minimised.pdb, a2b1
#load /media/rhys/data11/ampk_rhys_2022/metadynamics/data/a2b2+A769/01-Min/a2b2+A769_minimised.pdb, a2b2

#align a2b1 and resi 5-360 and name ca, a2b2 and resi 5-360 and name ca

# MAKE THE OLD FUNNEL
select s0_1, a2b1 and (resi 81 or resi 11) and name ca
select s1_1, a2b1 and (resi 302 or resi 43) and name ca
select s0_2, a2b2 and (resi 81 or resi 11) and name ca
select s1_2, a2b2 and (resi 303 or resi 43) and name ca

com s0_1, object=p0_1
com s1_1, object=p1_1
com s0_2, object=p0_2
com s1_2, object=p1_2

draw_funnel p0_1, p1=p0_1, p2=p1_1, s_cent=25, beta_cent=0.10, wall_width=15, wall_buffer=1.5, lower_wall=0, upper_wall=45, vec_step=2.5, angle_sample=18
set_name funnel, original_funnel_a2b1
draw_funnel p0_2, p1=p0_2, p2=p1_2, s_cent=25, beta_cent=0.10, wall_width=15, wall_buffer=1.5, lower_wall=0, upper_wall=45, vec_step=2.5, angle_sample=18
set_name funnel, original_funnel_a2b2


# MAKE THE NEW FUNNEL (OPEN SIMULATIONS)

select s0_1OP, a2b1 and (resi 81 or resi 11 or resi 310) and name ca
select s1_1OP, a2b1 and (resi 280 or resi 40) and name ca
select s0_2OP, a2b2 and (resi 81 or resi 11 or resi 311) and name ca
select s1_2OP, a2b2 and (resi 280 or resi 41) and name ca

com s0_1OP, object=p0_1OP
com s1_1OP, object=p1_1OP
com s0_2OP, object=p0_2OP
com s1_2OP, object=p1_2OP

draw_funnel p0_1OP, p1=p0_1OP, p2=p1_1OP, s_cent=25, beta_cent=0.15, wall_width=17.5, wall_buffer=2, lower_wall=-1, upper_wall=40, vec_step=2.5, angle_sample=18
set_name funnel, NEW_open_funnel_a2b1
draw_funnel p0_2OP, p1=p0_2OP, p2=p1_2OP, s_cent=25, beta_cent=0.15, wall_width=17.5, wall_buffer=2, lower_wall=-1, upper_wall=40, vec_step=2.5, angle_sample=18
set_name funnel, NEW_open_funnel_a2b2

draw_funnel p0_1OP, p1=p0_1OP, p2=p1_1OP, s_cent=25, beta_cent=0.15, wall_width=17.0, wall_buffer=2, lower_wall=-1, upper_wall=40, vec_step=2.5, angle_sample=18
set_name funnel, 17_NEW_open_funnel_a2b1
draw_funnel p0_2OP, p1=p0_2OP, p2=p1_2OP, s_cent=25, beta_cent=0.15, wall_width=17.0, wall_buffer=2, lower_wall=-1, upper_wall=40, vec_step=2.5, angle_sample=18
set_name funnel, 17_NEW_open_funnel_a2b2
