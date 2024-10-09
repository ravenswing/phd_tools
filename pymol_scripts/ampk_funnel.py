from pymol import cmd

systems = {
    "a2b1": ["A769", "MT47", "PF739", "MK87", "SC4"],
    "a2b2": ["PF739", "MK87", "SC4"],
}

pdb_dir = "/media/storage1/ampk_metad_all_data/open-fun-RMSD/structure_check"

cmd.set("orthoscopic", "on")

# MAKE THE NEW FUNNEL (OPEN SIMULATIONS)
for system, ligs in systems.items():
    for lig in ligs:
        cmd.load(f"{pdb_dir}/{system}+{lig}_R1_open.pdb", f"{system}+{lig}")

        cmd.select(f"{system}+{lig}_s0", f"{system}+{lig} and (resi 81 or resi 11 or resi 310) and name ca")
        cmd.select(f"{system}+{lig}_s1", f"{system}+{lig} and (resi 280 or resi 40) and name ca")

        cmd.do(f"com {system}+{lig}_s0, object={system}+{lig}_p0")
        cmd.do(f"com {system}+{lig}_s1, object={system}+{lig}_p1")

        cmd.delete(f"{system}+{lig}_s0")
        cmd.delete(f"{system}+{lig}_s1")

        cmd.do(f"draw_funnel {system}+{lig}_p0, p1={system}+{lig}_p0, p2={system}+{lig}_p1, s_cent=25, beta_cent=0.15, wall_width=17.5, wall_buffer=2, lower_wall=-1, upper_wall=40, vec_step=2.5, angle_sample=18")

        cmd.set_name( "funnel", f"{system}+{lig}_funnel")

cmd.show("spheres", "*_funnel")
cmd.set("sphere_scale", value=0.25, selection="*_funnel")

cmd.group("a2b1", "a2b1*")
cmd.group("a2b2", "a2b2*")
cmd.group("points", "a*_p*")



"""
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
"""
