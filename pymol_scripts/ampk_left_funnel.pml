

select s1_N1, a2b1 and (resi 301 or resi 302 or resi 13 or resi 12 or resi 25) and name ca
select s1_N2, a2b2 and (resi 301 or resi 302 or resi 13 or resi 12 or resi 25) and name ca

com s1_N1, object=p1_N1
com s1_N2, object=p1_N2

draw_funnel p0_1, p1_N1, s_cent=25, beta_cent=0.10, wall_width=15, wall_buffer=1.5, lower_wall=0, upper_wall=45, vec_step=2.5, angle_sample=18
set_name funnel, NEWfunnel_a2b1

draw_funnel p0_2, p1_N2, s_cent=25, beta_cent=0.10, wall_width=15, wall_buffer=1.5, lower_wall=0, upper_wall=45, vec_step=2.5, angle_sample=18
set_name funnel, NEWfunnel_a2b2

show spheres, NEWfunnel_a2b1
show spheres, NEWfunnel_a2b2

#hide lines, all
#show sticks, pdb1 and resn FRG
#show cartoon, pdb2
#show sticks, pdb2 and resn FRG

#color white, pdb1
#color lightorange, pdb2
#color magenta, (resn FRG)
#color forest, funnel
