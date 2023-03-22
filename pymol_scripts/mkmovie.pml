reinitialize
set movie_loop, 0

load movie.pdb

hide lines 
show sticks
hide (hydro)
spectrum count, rainbow, elem C
color gray70, (movie and resi 1-10 and elem C)
set_bond stick_transparency, 0.7, (movie and resi 11-120)

mset 1x1100

turn x, -90

frame 1
mview store
frame 30 
mview store


frame 80
turn x, 90
mview store


frame 99
mview store, state=1 

frame 1100
mview store, state=501 

mview reinterpolate, power=1

smooth
