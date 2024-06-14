from pymol import cmd
from glob import glob


def quick_vis(wd):
    for path in sorted(glob(f"{wd}/*.pdb")):
        name = path.split("/")[-1].split(".")[0]
        cmd.load(path, name)
        # cmd.show('cartoon')
        # cmd.rotate('x', -90, name)
        # cmd.rotate('y', -180, name)
