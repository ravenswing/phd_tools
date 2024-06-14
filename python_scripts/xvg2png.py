"""
===============================================================================
                        CONVERT XVG FILE TO PNG PLOT
===============================================================================
"""

import argparse
from glob import glob, iglob
import matplotlib.pyplot as plt
import subprocess


def convert_file(in_path, out_path):
    # Read in the xvg file
    with open(in_path) as f:
        lines = f.readlines()
    # Filter out comment lines
    data = [line for line in lines if line[0] not in ("@", "#")]
    # Convert data to numeric
    data = [[float(val) for val in line.split()[:2]] for line in data]
    # Seperate columns of data to plot
    x, y = [line[0] for line in data], [line[1] for line in data]
    # Create the simplest plot
    plt.plot(x, y)
    # Extract plot formatting info from the xvg output:
    for line in lines:
        # Title
        if line[0] == "@" and line.split()[1] == "title":
            plt.title(" ".join(line.split()[2:]).replace('"', ""))
        # X axis
        if line[0] == "@" and line.split()[1] == "xaxis":
            plt.xlabel(" ".join(line.split()[3:]).replace('"', ""))
        # Y axis
        if line[0] == "@" and line.split()[1] == "yaxis":
            plt.ylabel(" ".join(line.split()[3:]).replace('"', ""))
    # Save and close the new figure
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f" INFO | Successfully created: {pngname}")


if __name__ == "__main__":
    # Initiate argument parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=" "
    )
    # File or directory to process
    parser.add_argument("-f", type=str, help="File or Folder")
    # Optional arguments
    """
    parser.add_argument("-hisl", type=int, default=3,
                        help='Histidine in block (default: %(default)s)',
                        required=False)
    """
    # additional flags
    parser.add_argument(
        "-d",
        action="store_true",
        help=(
            "Display the plot once made - single file only." "(default: %(default)s)"
        ),
    )
    # Take variables from argparser
    args = parser.parse_args()
    PATH = args.f
    DISPLAY = args.d

    # Single file processing
    if ".xvg" in PATH:
        print(" INFO | OPERATING MODE - SINGLE FILE CONVERSION")
        print(f" INFO | Converting: {PATH}")
        # Define out_path
        pngname = f"{PATH[:-4]}.png"
        convert_file(PATH, pngname)

        # Display the file once it has been made
        if DISPLAY:
            try:
                subprocess.run(f"eog {pngname} &", check=True)
            except subprocess.CalledProcessError as error:
                print(
                    "Error code:",
                    error.returncode,
                    ". Output:",
                    error.output.decode("utf-8"),
                )
    # Whole directory processing
    else:
        print(" INFO | OPERATING MODE - WHOLE DIRECTORY CONVERSION")
        print(f" INFO | Converting: {PATH}")
        # Get the list of existing png files in that directory.
        pnglist = glob("./{}**/*.png".format(PATH), recursive=True)
        for filename in iglob("./{}**/*.xvg".format(PATH), recursive=True):
            pngname = f"{filename[:-4]}.png"
            # Skip any xvg files that have an existing png
            if pngname in pnglist:
                print(f" EXST | {pngname}")
                continue
            # Process the ones that do not and make new png.
            else:
                print(f" NEW  | {pngname}")
                convert_file(filename, pngname)
                plt.close()
