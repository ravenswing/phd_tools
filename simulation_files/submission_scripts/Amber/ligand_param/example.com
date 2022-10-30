%nproc=1
%chk=mol.chk
# b3lyp/6-31G(d) Opt

optimization mol

0 1

--Link1--
%nproc=1
%chk=mol.chk
# b3lyp/6-31G(d) POP=MK iop(6/50=1) geom=allcheck guess=read

mol.esp

