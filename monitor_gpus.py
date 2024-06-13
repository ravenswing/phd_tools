from fabric import Connection

# Make the connection to IQTC.
c = Connection(host='portal.qt.ub.edu', user='g19sllabres')

for n in [11, 12 ]:
    nvidia_out = c.run(f"ssh nodeg{n} -t \'nvidia-smi\'").stdout.split("\n")
    nvidia_out = 
    processes = c.run(f"ssh nodeg{n} -t \'ps -xh\'")
    print(len(processes.stdout.readlines()))
