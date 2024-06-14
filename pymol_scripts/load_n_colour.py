from pymol import cmd
from glob import glob


colours = {
    "a2b1-R1": "gray",
    "a2b1-R2": "brightorange",
    "a2b2-R1": "tv_blue",
    "a2b2-R2": "lightorange",
}


def quick_vis(wd):
    for path in sorted(glob(f"{wd}/*.pdb")):
        name = path.split("/")[-1].split(".")[0]
        c = name.split("+")[0] + "-" + name.split("_")[1]
        cmd.load(path, name)
        cmd.color(colours[c], selection=f"{name} and name C*")

        cmd.show("surface")
        cmd.refresh()


cmd.extend("load_clst", quick_vis)
