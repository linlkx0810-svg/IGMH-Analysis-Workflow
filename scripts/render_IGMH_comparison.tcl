if {[llength $argv] < 4} {
    puts "Usage: vmd -dispdev text -e scripts/render_IGMH_comparison.tcl -args SYSTEM CUBE_DIR OUTPUT_TGA TRANSFORM_TCL"
    quit
}

set system_name [lindex $argv 0]
set cubedir [lindex $argv 1]
set outfile [lindex $argv 2]
set transform_tcl [lindex $argv 3]

if {![file exists $transform_tcl]} {
    error "Comparison transform file not found: $transform_tcl"
}
source $transform_tcl

set img_w 2400
set img_h 2400
set sl2r_min -0.05
set sl2r_max  0.05
set dg_iso    0.005

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
light 2 off
light 3 off

color change rgb cyan 0.16 0.72 0.74
color Element H silver
color change rgb silver 0.78 0.78 0.76
color change rgb orange 0.72 0.42 0.06
color change rgb yellow 0.62 0.58 0.24
color Display Background white

material change ambient Opaque 0.14
material change diffuse Opaque 0.58
material change specular Opaque 0.08
material change shininess Opaque 0.18

material change opacity Transparent 0.60
material change ambient Transparent 0.16
material change diffuse Transparent 0.55
material change specular Transparent 0.02
material change shininess Transparent 0.10

proc apply_comparison_reps {} {
    global dg_iso sl2r_min sl2r_max
    mol delrep 0 top

    mol selection "all"
    mol representation CPK 0.500000 0.170000 32.000000 28.000000
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

proc apply_comparison_view {target} {
    global comparison_global comparison_center comparison_camera comparison_scale comparison_translate comparison_translate_panel comparison_projection
    if {![info exists comparison_global($target)]} {
        error "No molecule-level transform found for $target"
    }
    display resetview
    display projection $comparison_projection
    molinfo top set global_matrix [list $comparison_global($target)]
    molinfo top set center_matrix [list $comparison_center]
    molinfo top set rotate_matrix [list $comparison_camera]
    scale to $comparison_scale
    if {[info exists comparison_translate_panel($target)]} {
        set panel_translate $comparison_translate_panel($target)
    } else {
        set panel_translate $comparison_translate
    }
    puts "Panel screen translation for $target: $panel_translate"
    translate by [lindex $panel_translate 0] [lindex $panel_translate 1] [lindex $panel_translate 2]
}

puts "Rendering aligned comparison panel $system_name from $cubedir"
puts "Using transform file: $transform_tcl"
mol delete all
mol new "$cubedir/sl2r.cub" type cube waitfor all
mol addfile "$cubedir/dg_inter.cub" type cube waitfor all
apply_comparison_reps
apply_comparison_view $system_name
display update
render TachyonInternal $outfile
quit
