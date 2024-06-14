"""
USAGE:
ligand_res_names = {'A769': 'MOL',
                    'PF739':'PF',
                    'SC4':  'SC',
                    'MT47': 'MOL',
                    'MK87': 'MOL'}

# sys_info = list(product(*[['a2b1', 'a2b2'],['PF739'],['R1', 'R2', 'R3', 'R4']]))

for method in ['fun-metaD']:


    sys_info = list(product(*[['a2b1', 'a2b2'], ['A769', 'MT47', 'MK87'], ['R1', 'R2', 'R3', 'R4']]))

    tt.calculate_rmsd(sys_info,
                      f"{DATA_DIR}/{method}"+"/{p[0]}+{p[1]}/{p[2]}/md_dry.pdb",
                      f"{DATA_DIR}/{method}"+"/{p[0]}+{p[1]}/{p[2]}/metad_{p[0]}+{p[1]}_final.xtc",
                      f"{DATA_DIR}/analysis_data/{method}_ligand_rmsd.h5",
                      "resname MOL and not name H*",
                      ref_frmt=0,
                      unique_ligs=True)

    sys_info = list(product(*[['a2b1', 'a2b2'],['PF739'],['R1', 'R2', 'R3', 'R4']]))

    tt.calculate_rmsd(sys_info,
                      f"{DATA_DIR}/{method}"+"/{p[0]}+{p[1]}/{p[2]}/md_dry.pdb",
                      f"{DATA_DIR}/{method}"+"/{p[0]}+{p[1]}/{p[2]}/metad_{p[0]}+{p[1]}_final.xtc",
                      f"{DATA_DIR}/analysis_data/{method}_ligand_rmsd.h5",
                      "resname PF and not name H*",
                      ref_frmt=0,
                      unique_ligs=True)

    sys_info = list(product(*[['a2b1', 'a2b2'],['SC4'],['R1', 'R2', 'R3', 'R4']]))

    tt.calculate_rmsd(sys_info,
                      f"{DATA_DIR}/{method}"+"/{p[0]}+{p[1]}/{p[2]}/md_dry.pdb",
                      f"{DATA_DIR}/{method}"+"/{p[0]}+{p[1]}/{p[2]}/metad_{p[0]}+{p[1]}_final.xtc",
                      f"{DATA_DIR}/analysis_data/{method}_ligand_rmsd.h5",
                      "resname SC and not name H*",
                      ref_frmt=0,
                      unique_ligs=True)
"""


def calculate_rmsd(
    DIVS,
    top_frmt,
    trj_frmt,
    hdf_path,
    measure,
    ref_frmt=None,
    align="backbone",
    unique_ligs=False,
):
    all_sys = DIVS if unique_ligs else product(*DIVS)
    for i, p in enumerate(all_sys):
        p = list(p)

        print(f"Processing: {' '.join(p)}")

        if ref_frmt is not None:
            if isinstance(ref_frmt, int):
                ref = ref_frmt
            else:
                ref = ref_frmt.format(p=p)
        else:
            ref = top_frmt.format(p=p)

        print(f"Using Ref.: {ref} ({type(ref)})")

        # LIGAND & BACKBONE RMSD
        new_data = measure_rmsd(
            top_frmt.format(p=p), trj_frmt.format(p=p), ref, [measure], aln_group=align
        ).run()
        if len(all_sys[0]) == 3:
            inp = pd.DataFrame(
                columns=["t", p[2]], data=new_data.results.rmsd[:, [1, -1]]
            ).set_index("t")
            inp_l = pd.concat({p[1]: inp}, axis=1)
            inp_s = pd.concat({p[0]: inp_l}, axis=1)

        elif len(all_sys[0]) == 2:
            inp = pd.DataFrame(
                columns=["t", p[1]], data=new_data.results.rmsd[:, [1, -1]]
            ).set_index("t")
            inp_s = pd.concat({p[0]: inp}, axis=1)

        # If there is already an hdf file
        if exists(hdf_path):
            print("Further time --> Reading Files & Adding Data")
            new = pd.read_hdf(hdf_path, key="df")
            # Update the values if the data already exists
            if any([(mi == inp_s.columns)[0] for mi in new.columns]):
                print("Updating values in DataFrame.")
                new.update(inp_s)
            # Or add the new data
            else:
                print("Adding new values to DataFrame.")
                new = new.join(inp_s)
            # Reorder the columns before saving the data
            new = new.iloc[:, new.columns.sortlevel(0, sort_remaining=True)[1]]
            # Write the new data to the existing file
            new.to_hdf(hdf_path, key="df")

        # Or create one
        else:
            # Make a new hdf file and save the first column of data
            print("First time --> Creating Files")
            inp_s.to_hdf(hdf_path, key="df")
            i += 1
            continue


def calculate_rgyr(DIVS, top_frmt, trj_frmt, hdf_path, measure="protein"):
    for i, p in enumerate(product(*DIVS)):
        p = list(p)

        # LIGAND & BACKBONE RMSD
        new_data = measure_rgyr(top_frmt.format(p=p), trj_frmt.format(p=p), measure)
        if len(DIVS) == 3:
            inp = pd.DataFrame(columns=["t", p[2]], data=new_data).set_index("t")
            inp_l = pd.concat({p[1]: inp}, axis=1)
            inp_s = pd.concat({p[0]: inp_l}, axis=1)

        elif len(DIVS) == 2:
            inp = pd.DataFrame(
                columns=["t", p[1]], data=new_data.results.rmsd[:, [1, -1]]
            ).set_index("t")
            inp_s = pd.concat({p[0]: inp}, axis=1)

        # If there is already an hdf file
        if exists(hdf_path):
            print("Further time --> Reading Files & Adding Data")
            new = pd.read_hdf(hdf_path, key="df")
            # Update the values if the data already exists
            if any([(mi == inp_s.columns)[0] for mi in new.columns]):
                print("Updating values in DataFrame.")
                new.update(inp_s)
            # Or add the new data
            else:
                print("Adding new values to DataFrame.")
                new = new.join(inp_s)
            # Reorder the columns before saving the data
            new = new.iloc[:, new.columns.sortlevel(0, sort_remaining=True)[1]]
            # Write the new data to the existing file
            new.to_hdf(hdf_path, key="df")

        # Or create one
        else:
            # Make a new hdf file and save the first column of data
            print("First time --> Creating Files")
            inp_s.to_hdf(hdf_path, key="df")
            i += 1
            continue
