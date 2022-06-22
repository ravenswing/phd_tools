import argparse
import glob
import matplotlib.pyplot as plt
import subprocess

########################################################
#               GENERATE A-A SEQUENCES
########################################################

parser = argparse.ArgumentParser(
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        epilog=" ")

parser.add_argument("-f", type=str,
                    help='File or Folder')
# optional arguments
#parser.add_argument("-hisl", type=int, default=3,
                    #help='Histidine Lenght in block (default: %(default)s)',
                    #required=False)
# additional flags
parser.add_argument("-d", action="store_true",
                    help='Display once made single only(default: %(default)s)')

args = parser.parse_args()
PATH = args.f
DISPLAY = args.d

print(PATH)

if '.xvg' in PATH:
    print(" INFO | SINGLE FILE CONVERSION")
    print(f" INFO | Converting: {PATH}")
    pngname = f"{PATH[:-4]}.png"
    with open(PATH) as f:
        lines = f.readlines()
    data = [l for l in lines if l[0] not in ("@", "#")]
    data = [[float(val) for val in line.split()[:2]] for line in data]
    x, y = [l[0] for l in data], [l[1] for l in data]
    plt.plot(x, y)
    for line in lines:
        if line[0] == "@" and line.split()[1] == "title":
            plt.title(" ".join(line.split()[2:]).replace('"', ""))
        if line[0] == "@" and line.split()[1] == "xaxis":
            plt.xlabel(" ".join(line.split()[3:]).replace('"', ""))
        if line[0] == "@" and line.split()[1] == "yaxis":
            plt.ylabel(" ".join(line.split()[3:]).replace('"', ""))
    plt.savefig(pngname, dpi=300)
    print(f" INFO | Successfully created: {PATH[:-4]}.png")

    if DISPLAY:
        try:
            sub = subprocess.Popen(f"eog {pngname} &",
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT,
                                   shell=True
                                   )
            out, errors = sub.communicate()
        except:
            print(' ERROR | unable to display plot')

else:
    print(" INFO | WHOLE DIRECTORY CONVERSION")
    print(f" INFO | Converting: {PATH}")
    pnglist = glob.glob('./{}**/*.png'.format(PATH), recursive=True)
    for filename in glob.iglob('./{}**/*.xvg'.format(PATH), recursive=True):
        pngname = f"{filename[:-4]}.png"
        if pngname in pnglist:
            print(f" EXST | {pngname}")
            continue
        else:
            print(f" NEW  | {pngname}")
            with open(filename) as f:
                lines = f.readlines()
            data = [l for l in lines if l[0] not in ("@", "#")]
            data = [[float(val) for val in line.split()[:2]] for line in data]
            x, y = [l[0] for l in data], [l[1] for l in data]
            plt.plot(x, y)
            for line in lines:
                if line[0] == "@" and line.split()[1] == "title":
                    plt.title(" ".join(line.split()[2:]).replace('"', ""))
                if line[0] == "@" and line.split()[1] == "xaxis":
                    plt.xlabel(" ".join(line.split()[3:]).replace('"', ""))
                if line[0] == "@" and line.split()[1] == "yaxis":
                    plt.ylabel(" ".join(line.split()[3:]).replace('"', ""))
            plt.savefig(pngname)
            print(f" INFO | Successfully created: {pngname}")
            plt.close()
