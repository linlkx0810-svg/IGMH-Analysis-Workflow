if {[llength $argv] < 3} {
    puts "Usage: vmd -dispdev text -e scripts/render_IGMH.tcl -args SYSTEM CUBE_DIR OUTPUT_TGA"
    puts "Example: vmd -dispdev text -e scripts/render_IGMH.tcl -args TS3a output/TS3a_IGMH_files figures/TS3a_IGMH_clean.tga"
    quit
}

set system_name [lindex $argv 0]
set cubedir [lindex $argv 1]
set outfile [lindex $argv 2]
set center_selection ""
if {[llength $argv] > 3} {
    set center_selection [lindex $argv 3]
}

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

# Keep the IGMH blue/green/red scale intact; only soften molecular colors.
color change rgb cyan 0.16 0.72 0.74
color Element H silver
color change rgb silver 0.78 0.78 0.76
color change rgb orange 0.72 0.42 0.06
color change rgb yellow 0.62 0.58 0.24
color Display Background white

# Lower ambient/diffuse/specular terms for a calmer Tachyon render.
material change ambient Opaque 0.14
material change diffuse Opaque 0.58
material change specular Opaque 0.08
material change shininess Opaque 0.18

material change opacity Transparent 0.60
material change ambient Transparent 0.16
material change diffuse Transparent 0.55
material change specular Transparent 0.02
material change shininess Transparent 0.10

proc apply_refined_reps {} {
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

proc apply_refined_view {target} {
    display resetview
    rotate x by 62
    rotate y by -28
    rotate z by 34
    scale to 0.213

    if {$target eq "TS3a"} {
        translate by 0.05 -0.01 0.00
    } elseif {$target eq "TS3b"} {
        translate by 0.01 0.02 0.00
    }
}

proc center_on_selection {selection_text} {
    if {$selection_text eq ""} {
        return
    }
    set iface [atomselect top $selection_text]
    if {[$iface num] == 0} {
        puts "Warning: center_selection matched no atoms: $selection_text"
        $iface delete
        return
    }
    set c [measure center $iface weight none]
    $iface delete
    translate to 0 0 0
    translate by [vecscale -1.0 $c]
}

proc render_one {target cubedir outfile center_selection} {
    mol delete all
    mol new "$cubedir/sl2r.cub" type cube waitfor all
    mol addfile "$cubedir/dg_inter.cub" type cube waitfor all
    apply_refined_reps
    apply_refined_view $target
    center_on_selection $center_selection
    display update
    render TachyonInternal $outfile
}

puts "Rendering $system_name from $cubedir"
if {$center_selection ne ""} {
    puts "Center selection: $center_selection"
}
render_one $system_name $cubedir $outfile $center_selection
quit
