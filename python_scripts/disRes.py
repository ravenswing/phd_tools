import glob
import subprocess

LENGTH=16

gmx = '/usr/local/gromacs/bin/gmx'

stem = glob.glob('*.gro')[0].split('_')[0]
print(stem)

gro_file = stem+"_initial.gro"
index = "i.ndx"

try:
    subprocess.call("echo \"q\" | {G} make_ndx -f {gro} -o {ndx}".format(\
            G=gmx,gro=gro_file,ndx=index), \
            shell=True)
except:
    print('ERROR:make_ndx 1st run')

with open(index) as f:
    lines = f.readlines()
    ndx_end = len([l for l in lines if l[0] in ("[")])


for i in list(range(9)):
    make_ndx_command = "echo \" 4 & ri {}-{}\\n q\" | {G} make_ndx -f {gro} -n {ndx} -o {ndx}".format(\
            (i*LENGTH)+1, (i+1)*LENGTH, G=gmx,gro=gro_file,ndx=index)
    print(make_ndx_command)
    
    try:
        subprocess.call(make_ndx_command, shell=True)
    except:
        print('ERROR:make_ndx')
    itp_name = "disres_peptide{}.itp".format(i)
    genrestr_command = 'echo {} | {G} genrestr -f {gro} -n {ndx} -disre -o {itp}'.format(\
            i+ndx_end, G=gmx,gro=gro_file,ndx=index, itp=itp_name)
    try: 
        subprocess.call(genrestr_command, shell=True)
    except:
        print('ERROR:genrestr')

try:
    subprocess.call("cat disres_peptide0.itp > Disres.itp", shell=True)
except:
    print('ERROR:cat1')

for i in list(range(9))[1:]:
    try:
        subprocess.call("tail -n +5 disres_peptide{}.itp >> Disres.itp".format(i), shell=True)
    except:
        print('ERROR:cat2') 

try:
    subprocess.call("rm \\#*", shell=True)
except:
    print('ERROR:rm')

