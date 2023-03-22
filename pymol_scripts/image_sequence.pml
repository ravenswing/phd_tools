
reinitialize
run pymol_scripts/build_seq.py

build_seq AEAKAKAEAEAKKAEA, ss=helix
# build_seq EEEEKKKKEEEEKKKK, ss=helix
# build_seq EKEKEKEKEK, ss=helix


hide lines 
show cartoon
show sticks,(sidechain or name ca)
hide (hydro)

# set_color bright_yel, [ 255, 192, 0]
# color bright_yel, (name C*)

set_color royal_blu, [ 5, 4, 170]
color royal_blu, (name C*)

set ray_opaque_background, 0

# png /media/rhys/ExtHD/Project/carlos_peptides/LONG/hydrophilic_brush/helical_brush/figures/E4K4-1.png, width=1200, height=1200, dpi=300, ray=1

