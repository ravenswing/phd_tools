# Convert INPUT files from Amber -> Gromacs
import parmed as pmd


mols = [ 'H' ]
direc = 'Hydrophobic/'


for mol in mols:
    traj = pmd.load_file('{}{}/{}.prmtop'.format(direc,mol,mol), '{}{}/{}.rst7'.format(direc,mol,mol))
    traj.save('{}{}/{}.gro'.format(direc,mol,mol))
    traj.save('{}{}/{}.top'.format(direc,mol,mol))
