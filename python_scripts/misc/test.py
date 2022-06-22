import glob

paths = glob.glob('./Results/Hydrophobic/HHHG/*')
select = [ line[:-4].split('/') for line in paths if line[0] != '#' and '.edr' in line]

stem =  sorted(select,reverse=True)[0][-1]

print( stem )

print( stem[:-1]+str( int( stem[-1] ) +1) )

di = {'key':'value','key2':'value2'}

print( di[1] )
print( di[1][0])
