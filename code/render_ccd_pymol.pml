reinitialize
load data/5WN9_chainA_dnds_bfactor.pdb, ccd

bg_color white
set ray_opaque_background, 1
set orthoscopic, 1
set ray_shadows, 0
set antialias, 2
set specular, 0.25
set ambient, 0.45
set direct, 0.55
set cartoon_fancy_helices, 1
set label_size, 20
set label_color, black
set label_outline_color, white
set label_font_id, 7

hide everything
dss ccd

show cartoon, ccd
set cartoon_tube_radius, 0.45
cartoon tube, ccd
spectrum b, blue_white_red, ccd, minimum=-3, maximum=3

select noose, ccd and resi 173+176+182+186
show sticks, noose and sidechain
set stick_radius, 0.22, noose
color grey70, noose and elem C
color gold, noose and elem S
bond ccd and resi 173 and name SG, ccd and resi 186 and name SG
bond ccd and resi 176 and name SG, ccd and resi 182 and name SG
select ssbonds, (ccd and resi 173+186 and name SG) or (ccd and resi 176+182 and name SG)
show sticks, ssbonds
set stick_radius, 0.30, ssbonds

select cx3c, ccd and resi 182-186 and name CA
show spheres, cx3c
set sphere_scale, 0.30, cx3c

python
from pymol import cmd as _cmd
for _resi, _txt, _off in [("169", "N-169", (4, 4, 0)),
                          ("189", "C-189", (4, -4, 0)),
                          ("184", "CX3C 182-186", (-2, 8, 0))]:
    _x, _y, _z = _cmd.get_atom_coords("ccd and resi %s and name CA" % _resi)
    _cmd.pseudoatom("lbls", pos=[_x + _off[0], _y + _off[1], _z + _off[2]], label=_txt)
python end
hide nonbonded, lbls
set label_size, 24
set label_color, black
set label_outline_color, white

orient ccd
turn y, 20
turn x, -10
zoom (ccd or lbls), 2
ray 2400, 2000
png figures/rsv_g_ccd_structure_pymol_raw.png, dpi=300
