set research_dir "D:/Desktop/Research"
set ts3a_dir "$research_dir/TS3a_IGMH_files"
set ts3b_dir "$research_dir/TS3b_IGMH_files"

set img_w 2200
set img_h 2200
set sl2r_min -0.05
set sl2r_max  0.05
set dg_iso    0.005

# Fixed orthographic, publication-style display settings.
display resize $img_w $img_h
display projection Orthographic
display depthcue off
display shadows off
display ambientocclusion off
display antialias on
color Display Background white
axes location Off
color scale method BGR
light 0 on
light 1 on
light 2 on
light 3 on

# Thin, clean molecular representation; translucent IGMH surface.
material change ambient Opaque 0.22
material change diffuse Opaque 0.72
material change specular Opaque 0.18
material change shininess Opaque 0.45

material change opacity Transparent 0.60
material change ambient Transparent 0.28
material change diffuse Transparent 0.70
material change specular Transparent 0.08
material change shininess Transparent 0.25

proc apply_clean_reps {} {
    global dg_iso sl2r_min sl2r_max
    mol delrep 0 top

    mol selection "all"
    mol representation CPK 0.240000 0.070000 28.000000 22.000000
    mol color Element
    mol material Opaque
    mol addrep top

    mol selection "all"
    mol representation Isosurface $dg_iso 1 0 0 1 1
    mol color Volume 0
    mol material Transparent
    mol addrep top
    mol scaleminmax top 1 $sl2r_min $sl2r_max
}

proc center_on_interface {} {
    # VMD indices are zero-based. This centers the display on the interface:
    # fragment 2 atoms 73-94 plus the leading fragment-1 contributors near Co.
    set iface [atomselect top "index 0 1 2 11 70 71 72 to 93"]
    set c [measure center $iface weight none]
    $iface delete
    molinfo top set center_matrix [list [transoffset [vecscale -1.0 $c]]]
}

proc apply_clean_view {} {
    display resetview
    rotate x by 62
    rotate y by -28
    rotate z by 34
    scale to 0.195
}

proc render_one {cubedir outfile} {
    mol delete all
    mol new "$cubedir/sl2r.cub" type cube waitfor all
    mol addfile "$cubedir/dg_inter.cub" type cube waitfor all
    apply_clean_reps
    apply_clean_view
    display update
    render TachyonInternal $outfile
}

set target "both"
if {[llength $argv] > 0} {
    set target [lindex $argv 0]
}

if {$target eq "TS3a" || $target eq "both"} {
    render_one $ts3a_dir "$research_dir/TS3a_IGMH_clean.tga"
}
if {$target eq "TS3b" || $target eq "both"} {
    render_one $ts3b_dir "$research_dir/TS3b_IGMH_clean.tga"
}
quit
