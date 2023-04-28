delete all 

load /home/rhys/Dropbox/RESEARCH/AA_RHYS/BB_AMPK/Fun-metaD_Simulation_Files/System_Setup/apo/a2b1_apo.pdb, a2b1
load /home/rhys/Dropbox/RESEARCH/AA_RHYS/BB_AMPK/Fun-metaD_Simulation_Files/System_Setup/apo/a2b2_apo.pdb, a2b2

align a2b2 and resi 5-360 and name ca, a2b1 and resi 5-360 and name ca

select s0_1, a2b1 and (resi 81 and name ca or resi 11 and name ca) 
select s1_1, a2b1 and (resi 302 and name ca or resi 43 and name ca) 
select s0_2, a2b2 and (resi 81 and name ca or resi 11 and name ca) 
select s1_2, a2b2 and (resi 303 and name ca or resi 43 and name ca)

com s0_1, object=p0_1
com s1_1, object=p1_1
com s0_2, object=p0_2
com s1_2, object=p1_2

draw_funnel p0_1, p1_1, s_cent=25, beta_cent=0.10, wall_width=15, wall_buffer=1.5, lower_wall=0, upper_wall=45, vec_step=2.5, angle_sample=18
set_name funnel, funnel_a2b1
draw_funnel p0_2, p1_2, s_cent=25, beta_cent=0.10, wall_width=15, wall_buffer=1.5, lower_wall=0, upper_wall=45, vec_step=2.5, angle_sample=18
set_name funnel, funnel_a2b2

load /home/rhys/Storage/ampk_metad_all_data/a2b1+A769/0345-EQ-MD/a2b1+A769_lastframe.pdb, a2b1_A769
align a2b1_A769 and resi 5-360 and name ca, a2b1 and resi 5-360 and name ca
load /home/rhys/Storage/ampk_metad_all_data/a2b2+A769/0345-EQ-MD/a2b2+A769_lastframe.pdb, a2b2_A769
align a2b2_A769 and resi 5-360 and name ca, a2b1 and resi 5-360 and name ca

load /home/rhys/Storage/ampk_metad_all_data/a2b1+PF739/0345-EQ-MD/a2b1+PF739_lastframe.pdb, a2b1_PF739
align a2b1_PF739 and resi 5-360 and name ca, a2b1 and resi 5-360 and name ca
load /home/rhys/Storage/ampk_metad_all_data/a2b2+PF739/0345-EQ-MD/a2b2+PF739_lastframe.pdb, a2b2_PF739
align a2b2_PF739 and resi 5-360 and name ca, a2b1 and resi 5-360 and name ca

load /home/rhys/Storage/ampk_metad_all_data/a2b1+SC4/0345-EQ-MD/a2b1+SC4_lastframe.pdb, a2b1_SC4
align a2b1_SC4 and resi 5-360 and name ca, a2b1 and resi 5-360 and name ca
load /home/rhys/Storage/ampk_metad_all_data/a2b2+SC4/0345-EQ-MD/a2b2+SC4_lastframe.pdb, a2b2_SC4
align a2b2_SC4 and resi 5-360 and name ca, a2b1 and resi 5-360 and name ca

load /home/rhys/Storage/ampk_metad_all_data/a2b1+MT47/0345-EQ-MD/a2b1+MT47_lastframe.pdb, a2b1_MT47
align a2b1_MT47 and resi 5-360 and name ca, a2b1 and resi 5-360 and name ca
load /home/rhys/Storage/ampk_metad_all_data/a2b2+MT47/0345-EQ-MD/a2b2+MT47_lastframe.pdb, a2b2_MT47
align a2b2_MT47 and resi 5-360 and name ca, a2b1 and resi 5-360 and name ca

load /home/rhys/Storage/ampk_metad_all_data/a2b1+MK87/0345-EQ-MD/a2b1+MK87_lastframe.pdb, a2b1_MK87
align a2b1_MK87 and resi 5-360 and name ca, a2b1 and resi 5-360 and name ca
load /home/rhys/Storage/ampk_metad_all_data/a2b2+MK87/0345-EQ-MD/a2b2+MK87_lastframe.pdb, a2b2_MK87
align a2b2_MK87 and resi 5-360 and name ca, a2b1 and resi 5-360 and name ca

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


