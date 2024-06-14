import argparse
import sys
import numpy as np
import pandas as pd

# initialise the argument parser
PARSER = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="Reweighting Scripts",
    epilog=" ",
)

PARSER.add_argument("-a", action="store_true", help="Show all (default: %(default)s)")

PARSER.add_argument(
    "-p",
    type=str,
    default="5ake",
    help="Display PDB only(default: %(default)s)",
    required=False,
)


ARGS = PARSER.parse_args()

DATA_FILE = "./dG_database.h5"

database = pd.read_hdf(DATA_FILE, key="deltag")
if ARGS.a:
    print(database)
elif ARGS.p:
    pdb = ARGS.p
    print(database[pdb])
