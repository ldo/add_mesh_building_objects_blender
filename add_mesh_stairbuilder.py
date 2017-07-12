# Stairs and railing creator script for Blender
#
# Creates a straight-run staircase with railings and stringer
# All components are optional and can be turned on and off by setting e.g. makeTreads=True or makeTreads=False
# Current values assume 1 blender unit = 1 metre
#
# Stringer will rest on lower landing and hang from upper landing
# Railings start on the lowest step and end on the upper landing
#
# Note: I'm not sure how to use recalcNormals so not all normals points ouwards.
#       Perhaps someone else can contribute this.
#
#-----------------------------------------------------------
#
#   @todo:
#   - Join separate stringer objects and then clean up the mesh.
#   - Generate left/right posts/railings/retainers separately with
#       option to disable just the left/right.
#   - Add wall railing type as an option for left/right
#   - Add different rail styles (profiles).  Select with enum.
#   - Would like to add additional staircase types.
#       - Spiral staircase
#       - "L" staircase
#       - "T" staircase
#
# ##### BEGIN GPL LICENSE BLOCK #####
#
#  Stairbuilder is for quick stair generation.
#  Copyright (C) 2010  Nick van Adium
#  Copyright (C) 2011  Paul Marshall
#  Copyright (C) 2017 Lawrence D'Oliveiro
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ##### END GPL LICENSE BLOCK #####

import math
import enum
from copy import \
    copy
import bpy
from bpy.props import \
    BoolProperty, \
    EnumProperty, \
    IntProperty, \
    FloatProperty
from bpy_extras import \
    object_utils
import mathutils
from mathutils import \
    Matrix
from mathutils.geometry import \
    intersect_line_plane, \
    intersect_line_line

#+
# Useful stuff
#-

# deg = math.pi / 180 # degrees/radians conversion factor # not needed
circle = 2 * math.pi # circles/radians conversion factor

vec = lambda x, y, z : mathutils.Vector([x, y, z])
  # save some extra brackets

class EnumPropItems(enum.Enum) :
    "base class for enumerations that can be passed to Blender’s EnumProperty" \
    " to construct a menu of the enumeration values. Subclasses need only contain" \
    " one or more enumeration item definitions in the form\n" \
    "\n" \
    "    «name» = («title», «description»)\n" \
    "\n" \
    " in order for the all_items() method to return a tuple of tuples that can be" \
    " passed directly to EnumProperty."

    def __init__(self, label, description) :
        self._value_ = (label, description)
        self.label = label
        self.description = description
    #end __init__

    @classmethod
    def all_items(celf) :
        return \
            tuple((item.name, item.label, item.description) for item in celf.__members__.values())
    #end all_items

#end EnumPropItems

#+
# Building the parts
#-

class STAIRTYPE(EnumPropItems) :
    "overall types of stairs."
    FREESTANDING = ("Freestanding", "Generate a freestanding staircase.")
    HOUSED_OPEN = ("Housed-Open", "Generate a housed-open staircase.")
    BOX = ("Box", "Generate a box staircase.")
    CIRCULAR = ("Circular", "Generate a circular or spiral staircase.")
#end STAIRTYPE

class TREADTYPE(EnumPropItems) :
    "types of stair treads."
    CLASSIC = ("Classic", "Generate wooden style treads")
    BASIC_STEEL = ("Basic Steel", "Generate common steel style treads")
    BAR_1 = ("Bar 1", "Generate bar/slat steel treads")
    BAR_2 = ("Bar 2", "Generate bar-grating steel treads")
    BAR_3 = ("Bar 3", "Generate bar-support steel treads")
#end TREADTYPE

class STRINGERTYPE(EnumPropItems) :
    "types of stair stringers."
    CLASSIC = ("Classic", "Generate a classic style stringer")
    I_BEAM = ("I-Beam", "Generate a steel I-beam stringer")
    C_BEAM = ("C-Beam", "Generate a C-channel style stringer")
#end STRINGERTYPE

class MeshMaker :
    "object for creating meshes given the verts and faces."

    def __init__(self, rise, run, N) :
        # rise -- height of each tread
        # run -- depth of each tread
        # N -- number of treads
        self.stop = float(N) * vec(run, 0, rise)
        self.slope = rise / run
        # identical quads for all objects which are parallelpipeds (except stringers and treads)
        self.faces = \
            [
                [0, 1, 3, 2],
                [0, 1, 5, 4],
                [0, 2, 6, 4],
                [4, 5, 7, 6],
                [2, 3, 7, 6],
                [1, 3, 7, 5],
            ]
        self.made_objects = []
    #end __init__

    def make_mesh(self, verts, faces, name) :
        # MeshMaker, MeshMaker, make me a mesh...
        mesh = bpy.data.meshes.new(name)
        mesh.from_pydata(verts, [], faces)
        mesh.use_auto_smooth = True
        mesh.update()
        mesh_obj = bpy.data.objects.new(name = name, object_data = mesh)
        bpy.context.scene.objects.link(mesh_obj)
        mesh_obj.location = (0, 0, 0) # relative to root
        self.made_objects.append(mesh_obj)
    #end make_mesh

    def make_ppd_mesh(self, verts, name) :
        # special case of make_mesh for parallellipiped shapes
        assert len(verts) == 8
        self.make_mesh(verts, self.faces, name)
    #end make_ppd_mesh

#end MeshMaker

def posts(mm, rise, stair_run, post_depth, post_width, tread_width, nr_posts, rail_height, rail_thickness, rEnable, lEnable) :
    "generates posts for the stairs. These are the vertical elements holding up the railings."
    # TODO: STAIRTYPE.CIRCULAR

    x1 = vec(0, 0, rail_height - rail_thickness) #rail start
    x2 = mm.stop + vec(0, 0, rail_height - rail_thickness) #rail stop
    post_spacing = vec((x2[0] - x1[0]) / float(nr_posts + 1), 0, 0) #spacing between posts

    def intersect(i, d) :
        # finds intersection point, x, for rail and post
        x3 = x1 + i * post_spacing + vec(d, d, d)
        x4 = x3 + vec(0, 0, x2[-1])
        a = x2 - x1
        b = x4 - x3
        c = x3 - x1
        cr_ab = a.cross(b)
        mag_cr_ab = cr_ab * cr_ab
        return x1 + a * (c.cross(b).dot(cr_ab) / mag_cr_ab)
    #end intersect

#begin posts
    for i in range(nr_posts + 2) :
        coords = []
        #intersections with rail
        coords.append(intersect(i, 0.0))
        coords.append(intersect(i, post_depth))
        #intersections with tread
        coords.append(vec(
                x1[0] + i * post_spacing[0],
                0,
                int(coords[0][0] / stair_run) * rise
            ))
        coords.append(coords[2] + vec(post_depth, 0, 0))
        #inner face
        for j in range(4) :
            coords.append(coords[j] + vec(0, post_width, 0))
        #end for
        if rEnable :
            mm.make_ppd_mesh(coords, 'posts')
        #end if
        if lEnable :
            #make post on other side of steps as well
            for j in coords :
                j += vec(0, tread_width - post_width, 0)
            #end for
            mm.make_ppd_mesh(coords, 'posts')
        #end if
    #end for
#end posts

def railings(mm, rail_width, rail_thickness, rail_height, tread_toe, post_width, post_depth, tread_width, rEnable, lEnable) :
    "generates railings for the stairs. These go across the tops of the posts."
    # TODO: STAIRTYPE.CIRCULAR

    start = vec(0, 0, rail_height - rail_thickness) #rail start
    stop = mm.stop + vec(0, 0, rail_height - rail_thickness) #rail stop

    #determine offset to include railing toe
    offset = vec(tread_toe, 0, tread_toe * mm.slope)
    coords = []
    coords.append(start - offset)
    coords.append \
      (
            stop
        +
            offset
        +
            vec
              (
                post_depth,
                0,
                post_depth * mm.slope
              )
      )
    coords.append(start - offset + vec(0, rail_width, 0))
    coords.append \
      (
            stop
        +
            offset
        +
            vec
              (
                post_depth,
                rail_width,
                post_depth * mm.slope
              )
      )
    for j in range(4) :
        coords.append(coords[j] + vec(0, 0, rail_thickness))
    #end for
    #centre over posts
    for j in coords :
        j += vec(0, 0.5 * (- rail_width + post_width), 0)
    #end for
    if rEnable :
        mm.make_ppd_mesh(coords, 'rails')
    #end if
    if lEnable :
        #make rail on other side
        for j in coords :
            j += vec(0, tread_width - post_width, 0)
        #end for
        mm.make_ppd_mesh(coords, 'rails')
    #end if
#end railings

def retainers(mm, retainer_width, retainer_height, post_width, tread_width, rail_height, nr_retainers, rEnable, lEnable) :
    "generates retainers for the stairs. These are the additional pieces parallel" \
    " to, and below, the railings."

    retainer_spacing = rail_height / float(nr_retainers + 1)

    for i in range(nr_retainers) :
        coords = []
        offset = (i + 1) * vec(0, 0, retainer_spacing)
        coords.append(offset)
        coords.append(mm.stop + offset)
        coords.append(offset + vec(0, retainer_width, 0))
        coords.append(mm.stop + offset + vec(0, retainer_width, 0))
        for j in range(4) :
            coords.append(coords[j] + vec(0, 0, retainer_height))
        #end for
        #centre in posts
        for j in coords :
            j += vec(0, 0.5 * (post_width - retainer_width), 0)
        #end for
        if rEnable :
            mm.make_ppd_mesh(coords, 'retainers')
        #end if
        if lEnable :
            #make retainer on other side
            for j in coords :
                j += vec(0, tread_width - post_width, 0)
            #end for
            mm.make_ppd_mesh(coords, 'retainers')
        #end if
    #end for
#end retainers

def stringer(mm, stair_type, stringer_type, stair_rise, stair_run, w, stringer_height, nr_treads, tread_height, tread_width, tread_toe, tread_overhang, tw, stringer_flange_thickness, tp, stringer_intersects_ground,
    nr_stringers = 1, distributed_stringers = False, notMulti = True, sections_per_slice = 4) :
    "generates stringers for the stairs. These are the supports that go under" \
    " the stairs."

    if notMulti :
        stringer_width = w / 100
    else :
        stringer_width = (tread_width * (w / 100)) / nr_stringers
    #end if
    stringer_web_thickness = stringer_width * (tw / 100)
    stringer_flange_taper = 1 - tp / 100

    def freestanding_classic() :
        if distributed_stringers or nr_stringers == 1 :
            offset = tread_width / (nr_stringers + 1) - stringer_width / 2
        else :
            offset = 0
        #end if
        faces = \
            [
                [0, 1, 3, 2],
                [1, 5, 3],
                [3, 5, 4],
                [6, 7, 9, 8],
                [7, 11, 9],
                [9, 11, 10],
                [0, 2, 8, 6],
                [0, 1, 7, 6],
                [1, 5, 11, 7],
                [2, 3, 9, 8],
                [3, 4, 10, 9],
                [4, 5, 11, 10],
            ]
        for i in range(nr_stringers) :
            for j in range(nr_treads) :
                coords = []
                coords.append(vec(0, offset, - stair_rise))
                coords.append(vec(stair_run, offset, - stair_rise))
                coords.append(vec(0, offset, - tread_height))
                coords.append(vec(stair_run, offset, - tread_height))
                coords.append(vec(stair_run, offset, 0))
                coords.append(vec(stair_run * 2, offset, 0))
                for k in range(6) :
                    coords.append(coords[k] + vec(0, stringer_width, 0))
                #end for
                for k in coords :
                    k += j * vec(stair_run, 0, stair_rise)
                #end for
                mm.make_mesh(coords, faces, 'stringer')
            #end for
            if distributed_stringers or nr_stringers == 1 :
                offset += tread_width / (nr_stringers + 1)
            else :
                offset += (tread_width - stringer_width) / (nr_stringers - 1)
            #end if
        #end for
    #end freestanding_classic

    def housed_open_classic() :
        coords = []
        coords.append(vec(- tread_toe, - stringer_width, - stair_rise))
        coords.append(vec(tread_height / mm.slope, - stringer_width, - stair_rise))
        coords.append(vec(- tread_toe, - stringer_width, 0))
        coords.append(vec(
                nr_treads * stair_run,
                - stringer_width,
                (nr_treads - 1) * stair_rise - tread_height
            ))
        coords.append(vec(nr_treads * stair_run, - stringer_width, nr_treads * stair_rise))
        coords.append(vec(
                nr_treads * stair_run - tread_toe,
                - stringer_width,
                nr_treads * stair_rise
            ))
        for i in range(6) :
            coords.append(coords[i] + vec(0, stringer_width, 0))
        #end for
        faces = \
            [
                [0, 1, 7, 6],
                [1, 3, 9, 7],
                [3, 4, 10, 9],
                [4, 10, 11, 5],
                [5, 11, 8, 2],
                [2, 8, 6, 0],
                [0, 1, 2],
                [1, 2, 5, 3],
                [3, 4, 5],
                [6, 7, 8],
                [7, 8, 11, 9],
                [9, 10, 11],
            ]
        mm.make_mesh(coords, faces, 'stringer')
        for i in coords :
            i += vec(0, stringer_width + tread_width, 0)
        #end for
        mm.make_mesh(coords, faces, 'stringer')
    #end housed_open_classic

    def i_beam() :
        mid = stringer_width / 2
        web = stringer_web_thickness / 2
        # Bottom of the stringer:
        baseZ = - stair_rise - tread_height - stringer_height
        # Top of the strigner:
        topZ = - stair_rise - tread_height
        # Vertical taper amount:
        taper = stringer_flange_thickness * stringer_flange_taper
        if distributed_stringers or nr_stringers == 1 :
            offset = tread_width / (nr_stringers + 1) - mid
        else :
            offset = 0
        #end if
        # taper < 100%:
        if stringer_flange_taper > 0 :
            faces = \
                [
                    [0, 1, 17, 16],
                    [1, 2, 18, 17],
                    [2, 3, 19, 18],
                    [3, 4, 20, 19],
                    [4, 5, 21, 20],
                    [5, 6, 22, 21],
                    [6, 7, 23, 22],
                    [7, 8, 24, 23],
                    [8, 9, 25, 24],
                    [9, 10, 26, 25],
                    [10, 11, 27, 26],
                    [11, 12, 28, 27],
                    [12, 13, 29, 28],
                    [13, 14, 30, 29],
                    [14, 15, 31, 30],
                    [15, 0, 16, 31],
                    [0, 1, 2, 15],
                    [2, 11, 14, 15],
                    [11, 12, 13, 14],
                    [2, 3, 10, 11],
                    [3, 4, 5, 6],
                    [3, 6, 7, 10],
                    [7, 8, 9, 10],
                    [16, 17, 18, 31],
                    [18, 27, 30, 31],
                    [27, 28, 29, 30],
                    [18, 19, 26, 27],
                    [19, 20, 21, 22],
                    [19, 22, 23, 26],
                    [23, 24, 25, 26],
                ]
            for i in range(nr_stringers) :
                coords = []
                coords.append(vec(0, offset,                baseZ))
                coords.append(vec(0, offset,                baseZ + taper))
                coords.append(vec(0, offset + (mid - web),  baseZ + stringer_flange_thickness))
                coords.append(vec(0, offset + (mid - web),  topZ - stringer_flange_thickness))
                coords.append(vec(0, offset,                topZ - taper))
                coords.append(vec(0, offset,                topZ))
                coords.append(vec(0, offset + (mid - web),  topZ))
                coords.append(vec(0, offset + (mid + web),  topZ))
                coords.append(vec(0, offset + stringer_width,       topZ))
                coords.append(vec(0, offset + stringer_width,       topZ - taper))
                coords.append(vec(0, offset + (mid + web),  topZ - stringer_flange_thickness))
                coords.append(vec(0, offset + (mid + web),  baseZ + stringer_flange_thickness))
                coords.append(vec(0, offset + stringer_width,       baseZ + taper))
                coords.append(vec(0, offset + stringer_width,       baseZ))
                coords.append(vec(0, offset + (mid + web),  baseZ))
                coords.append(vec(0, offset + (mid - web),  baseZ))
                for j in range(16) :
                    coords.append(coords[j]+vec(stair_run * nr_treads, 0, stair_rise * nr_treads))
                #end for
                # If the bottom meets the ground:
                #   Bottom be flat with the xy plane, but shifted down.
                #   Either project onto the plane along a vector (hard) or use the built in
                #       interest found in mathutils.geometry (easy).  Using intersect:
                if stringer_intersects_ground :
                    for j in range(16) :
                        coords[j] = intersect_line_plane \
                          (
                            coords[j],
                            coords[j + 16],
                            vec(0, 0, topZ),
                            vec(0, 0, 1)
                          )
                    #end for
                #end if
                mm.make_mesh(coords, faces, 'stringer')
                if distributed_stringers or nr_stringers == 1 :
                    offset += tread_width / (nr_stringers + 1)
                else :
                    offset += (tread_width - stringer_width) / (nr_stringers - 1)
                #end if
            #end for
        # taper = 100%:
        else :
            faces = \
                [
                    [0, 1, 9, 8],
                    [1, 2, 10, 9],
                    [2, 3, 11, 10],
                    [3, 4, 12, 11],
                    [4, 5, 13, 12],
                    [5, 6, 14, 13],
                    [6, 7, 15, 14],
                    [7, 0, 8, 15],
                    [0, 1, 6, 7],
                    [1, 2, 5, 6],
                    [2, 3, 4, 5],
                    [8, 9, 14, 15],
                    [9, 10, 13, 14],
                    [10, 11, 12, 13],
                ]
            for i in range(nr_stringers) :
                coords = []
                coords.append(vec(0, offset,                baseZ))
                coords.append(vec(0, offset + (mid - web),  baseZ + stringer_flange_thickness))
                coords.append(vec(0, offset + (mid - web),  topZ - stringer_flange_thickness))
                coords.append(vec(0, offset,                topZ))
                coords.append(vec(0, offset + stringer_width,       topZ))
                coords.append(vec(0, offset + (mid + web),  topZ - stringer_flange_thickness))
                coords.append(vec(0, offset + (mid + web),  baseZ + stringer_flange_thickness))
                coords.append(vec(0, offset + stringer_width,       baseZ))
                for j in range(8) :
                    coords.append(coords[j] + vec(stair_run * nr_treads, 0, stair_rise * nr_treads))
                #end for
                mm.make_mesh(coords, faces, 'stringer')
                offset += tread_width / (nr_stringers + 1)
            #end for
        #end if
    #end i_beam

    def housed_i_beam() :
        webOrth = vec(stair_rise, 0, - stair_run).normalized()
        webHeight = vec(stair_run + tread_toe, 0, - tread_height).project(webOrth).length
        vDelta_1 = stringer_flange_thickness * mm.slope
        vDelta_2 = stair_rise * (nr_treads - 1) - (webHeight + stringer_flange_thickness)
        flange_y = (stringer_width - stringer_web_thickness) / 2
        front = - tread_toe - stringer_flange_thickness
        outer = - tread_overhang - stringer_web_thickness - flange_y
        coords = []
        if stringer_flange_taper > 0 :
            # Upper-Outer flange:
            coords.append(vec(front, outer, - stair_rise))
            coords.append(vec(- tread_toe, outer, - stair_rise))
            coords.append(vec(- tread_toe, outer, 0))
            coords.append(vec(
                    stair_run * (nr_treads - 1) - tread_toe,
                    outer,
                    stair_rise * (nr_treads - 1)
                ))
            coords.append(vec(
                    stair_run * nr_treads,
                    outer,
                    stair_rise * (nr_treads - 1)
                ))
            coords.append(vec(
                    stair_run * nr_treads,
                    outer,
                    stair_rise * (nr_treads - 1) + stringer_flange_thickness
                ))
            coords.append(vec(
                    stair_run * (nr_treads - 1) - tread_toe,
                    outer,
                    stair_rise * (nr_treads - 1) + stringer_flange_thickness
                ))
            coords.append(vec(front, outer, stringer_flange_thickness - vDelta_1))
            # Lower-Outer flange:
            coords.append(coords[0] + vec(stringer_flange_thickness + webHeight, 0, 0))
            coords.append(coords[1] + vec(stringer_flange_thickness + webHeight, 0, 0))
            coords.append \
              (
                intersect_line_line
                  (
                    coords[9],
                    coords[9] - vec(0, 0, 1),
                    vec(stair_run, 0, - tread_height - stringer_flange_thickness),
                    vec(stair_run * 2, 0, stair_rise - tread_height - stringer_flange_thickness)
                 )[0]
              )
            coords.append(vec(
                    stair_run * nr_treads - (webHeight - tread_height) / mm.slope,
                    outer,
                    vDelta_2
                ))
            coords.append(coords[4] - vec(0, 0, stringer_flange_thickness + webHeight))
            coords.append(coords[5] - vec(0, 0, stringer_flange_thickness + webHeight))
            coords.append(coords[11] + vec(0, 0, stringer_flange_thickness))
            coords.append \
              (
                  intersect_line_line
                    (
                      coords[8],
                      coords[8] - vec(0, 0, 1),
                      vec(stair_run, 0, - tread_height),
                      vec(stair_run * 2, 0, stair_rise - tread_height)
                    )[0]
              )
            # Outer web:
            coords.append(coords[1] + vec(0, flange_y, 0))
            coords.append(coords[8] + vec(0, flange_y, 0))
            coords.append(coords[15] + vec(0, flange_y, 0))
            coords.append(coords[14] + vec(0, flange_y, 0))
            coords.append(coords[13] + vec(0, flange_y, 0))
            coords.append(coords[4] + vec(0, flange_y, 0))
            coords.append(coords[3] + vec(0, flange_y, 0))
            coords.append(coords[2] + vec(0, flange_y, 0))
            # Upper-Inner flange and lower-inner flange:
            for i in range(16) :
                coords.append(coords[i] + vec(0, stringer_width, 0))
            #end for
            # Inner web:
            for i in range(8) :
                coords.append(coords[i + 16] + vec(0, stringer_web_thickness, 0))
            #end for
            # Mid nodes to so faces will be quads:
            for i in [0, 7, 6, 5, 9, 10, 11, 12] :
                coords.append(coords[i] + vec(0, flange_y, 0))
            #end for
            for i in range(8) :
                coords.append(coords[i + 48] + vec(0, stringer_web_thickness, 0))
            #end for
            faces = \
                [
                    [0, 1, 2, 7],
                    [2, 3, 6, 7],
                    [3, 4, 5, 6],
                    [1, 2, 23, 16],
                    [2, 3, 22, 23],
                    [3, 4, 21, 22],
                    [16, 17, 18, 23],
                    [18, 19, 22, 23],
                    [19, 20, 21, 22],
                    [17, 8, 15, 18],
                    [18, 15, 14, 19],
                    [19, 14, 13, 20],
                    [8, 9, 10, 15],
                    [10, 11, 14, 15],
                    [11, 12, 13, 14],
                    [9, 10, 53, 52],
                    [10, 11, 54, 53],
                    [11, 12, 55, 54],
                    [52, 53, 61, 60],
                    [53, 54, 62, 61],
                    [54, 55, 63, 62],
                    [60, 61, 34, 33],
                    [61, 62, 35, 34],
                    [62, 63, 36, 35],
                    [32, 33, 34, 39],
                    [34, 35, 38, 39],
                    [35, 36, 37, 38],
                    [41, 32, 39, 42],
                    [42, 39, 38, 43],
                    [43, 38, 37, 44],
                    [40, 41, 42, 47],
                    [42, 43, 46, 47],
                    [43, 44, 45, 46],
                    [25, 26, 47, 40],
                    [26, 27, 46, 47],
                    [27, 28, 45, 46],
                    [24, 25, 26, 31],
                    [26, 27, 30, 31],
                    [27, 28, 29, 30],
                    [24, 31, 57, 56],
                    [31, 30, 58, 57],
                    [30, 29, 59, 58],
                    [48, 49, 57, 56],
                    [49, 50, 58, 57],
                    [50, 51, 59, 58],
                    [0, 7, 49, 48],
                    [7, 6, 50, 49],
                    [6, 5, 51, 50],
                    [0, 1, 16, 48],
                    [16, 40, 56, 48],
                    [24, 25, 40, 56],
                    [16, 17, 41, 40],
                    [8, 9, 52, 17],
                    [17, 52, 60, 41],
                    [32, 33, 60, 41],
                    [12, 13, 20, 55],
                    [20, 44, 63, 55],
                    [37, 44, 63, 36],
                    [20, 21, 45, 44],
                    [28, 29, 51, 21],
                    [21, 51, 59, 45],
                    [28, 45, 59, 29],
                    [4, 5, 51, 21],
                ]
            mm.make_mesh(coords, faces, 'stringer')
            for i in coords :
                i += vec(0, tread_width + stringer_web_thickness, 0)
            #end for
            mm.make_mesh(coords, faces, 'stringer')
        #end if
        # @TODO Taper = 100%
    #end housed_i_beam

    def housed_c_beam() :
        webOrth = vec(stair_rise, 0, - stair_run).normalized()
        webHeight = vec(stair_run + tread_toe, 0, - tread_height).project(webOrth).length
        vDelta_1 = stringer_flange_thickness * mm.slope
        vDelta_2 = stair_rise * (nr_treads - 1) - (webHeight + stringer_flange_thickness)
        flange_y = (stringer_width - stringer_web_thickness) / 2
        front = - tread_toe - stringer_flange_thickness
        outer = - tread_overhang - stringer_web_thickness - flange_y
        coords = []
        if stringer_flange_taper > 0 :
            # Upper-Outer flange:
            coords.append(vec(front, outer, - stair_rise))
            coords.append(vec(- tread_toe, outer, - stair_rise))
            coords.append(vec(- tread_toe, outer, 0))
            coords.append(vec(
                    stair_run * (nr_treads - 1) - tread_toe,
                    outer,
                    stair_rise * (nr_treads - 1)
                ))
            coords.append(vec(
                    stair_run * nr_treads,
                    outer,
                    stair_rise * (nr_treads - 1)
                ))
            coords.append(vec(
                    stair_run * nr_treads,
                    outer,
                    stair_rise * (nr_treads - 1) + stringer_flange_thickness
                ))
            coords.append(vec(
                    stair_run * (nr_treads - 1) - tread_toe,
                    outer,
                    stair_rise * (nr_treads - 1) + stringer_flange_thickness
                ))
            coords.append(vec(front, outer, stringer_flange_thickness - vDelta_1))
            # Lower-Outer flange:
            coords.append(coords[0] + vec(stringer_flange_thickness + webHeight, 0, 0))
            coords.append(coords[1] + vec(stringer_flange_thickness + webHeight, 0, 0))
            coords.append \
              (
                intersect_line_line
                  (
                    coords[9],
                    coords[9] - vec(0, 0, 1),
                    vec(stair_run, 0, - tread_height - stringer_flange_thickness),
                    vec(stair_run * 2, 0, stair_rise - tread_height - stringer_flange_thickness)
                )[0]
              )
            coords.append(vec(
                    stair_run * nr_treads - (webHeight - tread_height) / mm.slope,
                    outer,
                    vDelta_2
                ))
            coords.append(coords[4] - vec(0, 0, stringer_flange_thickness + webHeight))
            coords.append(coords[5] - vec(0, 0, stringer_flange_thickness + webHeight))
            coords.append(coords[11] + vec(0, 0, stringer_flange_thickness))
            coords.append \
              (
                intersect_line_line
                  (
                    coords[8],
                    coords[8] - vec(0, 0, 1),
                    vec(stair_run, 0, - tread_height),
                    vec(stair_run * 2, 0, stair_rise - tread_height)
                  )[0]
              )
            # Outer web:
            coords.append(coords[1] + vec(0, flange_y, 0))
            coords.append(coords[8] + vec(0, flange_y, 0))
            coords.append(coords[15] + vec(0, flange_y, 0))
            coords.append(coords[14] + vec(0, flange_y, 0))
            coords.append(coords[13] + vec(0, flange_y, 0))
            coords.append(coords[4] + vec(0, flange_y, 0))
            coords.append(coords[3] + vec(0, flange_y, 0))
            coords.append(coords[2] + vec(0, flange_y, 0))
            # Outer corner nodes:
            for i in [0, 7, 6, 5, 12, 11, 10, 9] :
                coords.append(coords[i] + vec(0, flange_y + stringer_web_thickness, 0))
            #end for
            faces = \
                [
                    [0, 1, 2, 7],
                    [2, 3, 6, 7],
                    [3, 4, 5, 6],
                    [1, 2, 23, 16],
                    [2, 3, 22, 23],
                    [3, 4, 21, 22],
                    [16, 17, 18, 23],
                    [18, 19, 22, 23],
                    [19, 20, 21, 22],
                    [17, 8, 15, 18],
                    [18, 15, 14, 19],
                    [19, 14, 13, 20],
                    [8, 9, 10, 15],
                    [10, 11, 14, 15],
                    [11, 12, 13, 14],
                    [0, 24, 25, 7],
                    [7, 25, 26, 6],
                    [6, 26, 27, 5],
                    [9, 31, 30, 10],
                    [10, 30, 29, 11],
                    [11, 29, 28, 12],
                    [24, 25, 30, 31],
                    [25, 26, 29, 30],
                    [26, 27, 28, 29],
                    [0, 1, 16, 24],
                    [16, 24, 31, 17],
                    [8, 9, 31, 17],
                    [4, 5, 27, 21],
                    [20, 21, 27, 28],
                    [12, 13, 20, 28],
                ]
            mm.make_mesh(coords, faces, 'stringer')
            for i in range(16) :
                coords[i] += vec(0, - outer * 2, 0)
            #end for
            for i in range(8) :
                coords[i + 16] += vec(0, (- outer - flange_y) * 2, 0)
            #end for
            for i in coords :
                i += vec(0, tread_overhang * 2 + tread_width, 0)
            #end for
            mm.make_mesh(coords, faces, 'stringer')
        #end if
    #end housed_c_beam

    def box() :
        h = - tread_height #height of top section
        for i in range(nr_treads) :
            coords = []
            coords.append(vec(i * stair_run, 0, - stair_rise))
            coords.append(vec((i + 1) * stair_run, 0, - stair_rise))
            coords.append(vec(i * stair_run, 0, h + i * stair_rise))
            coords.append(vec((i + 1) * stair_run, 0, h + i * stair_rise))
            for j in range(4) :
                coords.append(coords[j] + vec(0, tread_width, 0))
            #end for
            mm.make_ppd_mesh(coords, 'stringer')
        #end for
    #end box

    def circular() :
        offset = tread_width / (nr_stringers + 1) - stringer_width / 2
        for s in range(nr_stringers) :
            base = tread_overhang + offset * (s + 1)
            start = \
                [
                    vec(0, - base, - tread_height),
                    vec(0, - base, - tread_height - stair_rise),
                    vec(0, - base - stringer_width, - tread_height),
                    vec(0, - base - stringer_width, - tread_height - stair_rise),
                ]
            tread_angle = stair_run / nr_treads
            for i in range(nr_treads) :
                coords = []
                # Base faces.  Should be able to append more sections :
                bar_2_faces = [[0, 1, 3, 2]]
                t_inner = Matrix.Rotation(tread_angle * i, 3, 'Z')
                coords.append(t_inner * start[0] + vec(0, 0, stair_rise * i))
                coords.append(t_inner * start[1] + vec(0, 0, stair_rise * i))
                t_outer = Matrix.Rotation(tread_angle * i, 3, 'Z')
                coords.append(t_outer * start[2] + vec(0, 0, stair_rise * i))
                coords.append(t_outer * start[3] + vec(0, 0, stair_rise * i))
                for j in range(sections_per_slice) :
                    k = j * 4 + 4
                    bar_2_faces.append([k, k - 4, k - 3, k + 1])
                    bar_2_faces.append([k - 2, k - 1, k + 3, k + 2])
                    bar_2_faces.append([k + 1, k - 3, k - 1, k + 3])
                    bar_2_faces.append([k, k - 4, k - 2, k + 2])
                    rot = Matrix.Rotation \
                      (
                        tread_angle * (j + 1) / sections_per_slice + tread_angle * i,
                        3,
                        'Z'
                      )
                    for v in start :
                        coords.append(rot * v + vec(0, 0, stair_rise * i))
                    #end for
                #end for
                for j in range(sections_per_slice) :
                    k = (j + sections_per_slice) * 4 + 4
                    bar_2_faces.append([k, k - 4, k - 3, k + 1])
                    bar_2_faces.append([k - 2, k - 1, k + 3, k + 2])
                    bar_2_faces.append([k + 1, k - 3, k - 1, k + 3])
                    bar_2_faces.append([k, k - 4, k - 2, k + 2])
                    rot = Matrix.Rotation \
                      (
                        tread_angle * (j + sections_per_slice + 1) / sections_per_slice + tread_angle * i,
                        3,
                        'Z'
                      )
                    for v in range(4) :
                        if v in [1, 3] :
                            incline = stair_rise * i + stair_rise / sections_per_slice * (j + 1)
                            coords.append(rot * start[v] + vec(0, 0, incline))
                        else :
                            coords.append(rot * start[v] + vec(0, 0, stair_rise * i))
                        #end if
                    #end for
                #end for
                mm.make_mesh(coords, bar_2_faces, 'stringer')
            #end for
        #end for
    #end circular

#begin stringer
    if stair_type == STAIRTYPE.FREESTANDING :
        if stringer_type == STRINGERTYPE.CLASSIC :
            freestanding_classic()
        elif stringer_type == STRINGERTYPE.I_BEAM :
            i_beam()
        # note STRINGERTYPE.C_BEAM not supported
        #end if
    elif stair_type == STAIRTYPE.HOUSED_OPEN :
        if stringer_type == STRINGERTYPE.CLASSIC :
            housed_open_classic()
        elif stringer_type == STRINGERTYPE.I_BEAM :
            housed_i_beam()
        elif stringer_type == STRINGERTYPE.C_BEAM :
            housed_c_beam()
        #end if
    elif stair_type == STAIRTYPE.BOX :
        box()
    elif stair_type == STAIRTYPE.CIRCULAR :
        circular()
    #end if
#end stringer

def treads(mm, stair_type, tread_type, stair_run, tread_width, tread_height, tread_run, tread_rise, tread_toe, tread_side_overhang, nr_treads, tread_metal_thickness, nr_tread_sections, spacing, nr_cross_sections) :
    "generates treads for non-circular stairs."

    if nr_tread_sections != 1 and tread_type not in [TREADTYPE.BAR_2, TREADTYPE.BAR_3] :
        tread_section_spacing = (tread_run + tread_toe) * (spacing / 100) / (nr_tread_sections - 1)
          # spacing between sections (% of depth)
    elif tread_type in [TREADTYPE.BAR_2, TREADTYPE.BAR_3] :
        tread_section_spacing = spacing / 100 # keep % value
    else :
        tread_section_spacing = 0
    #end if

    # Setup the coordinates:
    treads_coords = []
    bars_coords = []
    crosses_coords = []
    cross = 0
    cW = 0
    depth = 0
    tread_section_offset = 0
    height = 0
    assert stair_type in [STAIRTYPE.FREESTANDING, STAIRTYPE.HOUSED_OPEN, STAIRTYPE.BOX]

    def make_treads_coords() :
        # calculates coordinates for the pieces of the treads.
        nonlocal tread_section_spacing, tread_section_offset, cross, cW, depth, height
        if tread_type == TREADTYPE.CLASSIC :
            treads_coords.append(vec(- tread_toe, - tread_side_overhang, 0))
            treads_coords.append(vec(tread_run, - tread_side_overhang, 0))
            treads_coords.append(vec(- tread_toe, tread_width + tread_side_overhang, 0))
            treads_coords.append(vec(tread_run, tread_width + tread_side_overhang, 0))
            for i in range(4) :
                treads_coords.append(treads_coords[i] + vec(0, 0, - tread_height))
            #end for
        elif tread_type == TREADTYPE.BASIC_STEEL :
            depth = (tread_run + tread_toe - (nr_tread_sections - 1) * tread_section_spacing) / nr_tread_sections
            inset = depth / 4
            tDepth = depth - tread_toe
            treads_coords.append(vec(- tread_toe, - tread_side_overhang, - tread_height))                          #0
            treads_coords.append(vec(inset - tread_toe, - tread_side_overhang, - tread_height))           #1
            treads_coords.append(vec(inset - tread_toe, - tread_side_overhang, - tread_height + tread_metal_thickness)) #2
            treads_coords.append(vec(tread_metal_thickness - tread_toe, - tread_side_overhang, - tread_height + tread_metal_thickness))       #3
            treads_coords.append(vec(tread_metal_thickness - tread_toe, - tread_side_overhang, - tread_metal_thickness))                #4
            treads_coords.append(vec(- tread_toe, - tread_side_overhang, 0))                                #5
            treads_coords.append(vec(tDepth, - tread_side_overhang, 0))                                 #6
            treads_coords.append(vec(tDepth - tread_metal_thickness, - tread_side_overhang, - tread_metal_thickness))                #7
            treads_coords.append(vec(tDepth - tread_metal_thickness, - tread_side_overhang, tread_metal_thickness - tread_height))        #8
            treads_coords.append(vec(tDepth, - tread_side_overhang, - tread_height))                           #9
            treads_coords.append(vec(tDepth - inset, - tread_side_overhang, - tread_height))           #10
            treads_coords.append(vec(tDepth - inset, - tread_side_overhang, - tread_height + tread_metal_thickness)) #11
            for i in range(12) :
                treads_coords.append(treads_coords[i] + vec(0, tread_width + 2 * tread_side_overhang, 0))
            #end for
        elif tread_type in [TREADTYPE.BAR_1, TREADTYPE.BAR_2, TREADTYPE.BAR_3] :
            # Frame:
            treads_coords.append(vec(- tread_toe, - tread_side_overhang, - tread_height))
            treads_coords.append(vec(tread_run, - tread_side_overhang, - tread_height))
            treads_coords.append(vec(- tread_toe, - tread_side_overhang, 0))
            treads_coords.append(vec(tread_run, - tread_side_overhang, 0))
            for i in range(4) :
                if i % 2 == 0 :
                    treads_coords.append(treads_coords[i] + vec(tread_metal_thickness, tread_metal_thickness, 0))
                else :
                    treads_coords.append(treads_coords[i] + vec(- tread_metal_thickness, tread_metal_thickness, 0))
                #end if
            #end for
            for i in range(4) :
                treads_coords.append(treads_coords[i] + vec(0, tread_width + tread_side_overhang, 0))
            #end for
            for i in range(4) :
                treads_coords.append(treads_coords[i + 4] + vec(0, tread_width + tread_side_overhang - 2 * tread_metal_thickness, 0))
            #end for
            # Tread sections:
            if tread_type == TREADTYPE.BAR_1 :
                tread_section_offset = (tread_metal_thickness * math.sqrt(2)) / 2
                topset = tread_height - tread_section_offset
                tread_section_spacing = ((tread_run + tread_toe - 2 * tread_metal_thickness) - (tread_section_offset * nr_tread_sections + topset)) / (nr_tread_sections + 1)
                baseX = - tread_toe + tread_section_spacing + tread_metal_thickness
                bars_coords.append(vec(baseX, tread_metal_thickness - tread_side_overhang, tread_section_offset - tread_height))
                bars_coords.append(vec(baseX + tread_section_offset, tread_metal_thickness - tread_side_overhang, - tread_height))
                for i in range(2) :
                    bars_coords.append(bars_coords[i] + vec(topset, 0, topset))
                #end for
                for i in range(4) :
                    bars_coords.append(bars_coords[i] + vec(0, (tread_width + tread_side_overhang) - 2 * tread_metal_thickness, 0))
                #end for
            elif tread_type in [TREADTYPE.BAR_2, TREADTYPE.BAR_3] :
                tread_section_offset = ((stair_run + tread_toe) * tread_section_spacing) / (nr_tread_sections + 1)
                topset = ((stair_run + tread_toe) * (1 - tread_section_spacing) - 2 * tread_metal_thickness) / nr_tread_sections
                baseX = - tread_toe + tread_metal_thickness + tread_section_offset
                baseY = tread_width + tread_side_overhang - 2 * tread_metal_thickness
                bars_coords.append(vec(baseX, - tread_side_overhang + tread_metal_thickness, - tread_height / 2))
                bars_coords.append(vec(baseX + topset, - tread_side_overhang + tread_metal_thickness, - tread_height / 2))
                bars_coords.append(vec(baseX, - tread_side_overhang + tread_metal_thickness, 0))
                bars_coords.append(vec(baseX + topset, - tread_side_overhang + tread_metal_thickness, 0))
                for i in range(4) :
                    bars_coords.append(bars_coords[i] + vec(0, baseY, 0))
                #end for
            #end if
            # Tread cross-sections:
            if tread_type in [TREADTYPE.BAR_1, TREADTYPE.BAR_2] :
                cW = tread_metal_thickness
                cross = (tread_width + 2 * tread_side_overhang - (nr_cross_sections + 2) * tread_metal_thickness) / (nr_cross_sections + 1)
            else : # TREADTYPE.BAR_3
                spacing = tread_section_spacing ** (1 / 4)
                cross = ((2 * tread_side_overhang + tread_width) * spacing) / (nr_cross_sections + 1)
                cW = (- 2 * tread_metal_thickness + (2 * tread_side_overhang + tread_width) * (1 - spacing)) / nr_cross_sections
                tread_section_spacing = topset
                height = - tread_height / 2
            #end if
            baseY = - tread_side_overhang + tread_metal_thickness + cross
            crosses_coords.append(vec(- tread_toe + tread_metal_thickness, baseY, - tread_height))
            crosses_coords.append(vec(tread_run - tread_metal_thickness, baseY, - tread_height))
            crosses_coords.append(vec(- tread_toe + tread_metal_thickness, baseY, height))
            crosses_coords.append(vec(tread_run - tread_metal_thickness, baseY, height))
            for i in range(4) :
                crosses_coords.append(crosses_coords[i] + vec(0, cW, 0))
            #end for
        #end if
    #end make_treads_coords

    def make_treads() :
        # actually creates the objects for the tread components.
        basic_steel_faces = \
            [
                [0, 1, 2, 3],
                [0, 3, 4, 5],
                [4, 5, 6, 7],
                [6, 7, 8, 9],
                [8, 9, 10, 11],
                [12, 13, 14, 15],
                [12, 15, 16, 17],
                [16, 17, 18, 19],
                [18, 19, 20, 21],
                [20, 21, 22, 23],
                [0, 1, 13, 12],
                [1, 2, 14, 13],
                [2, 3, 15, 14],
                [3, 4, 16, 15],
                [4, 7, 19, 16],
                [7, 8, 20, 19],
                [8, 11, 23, 20],
                [11, 10, 22, 23],
                [10, 9, 21, 22],
                [9, 6, 18, 21],
                [6, 5, 17, 18],
                [5, 0, 12, 17],
            ]
        bar_faces = \
            [
                [0 , 2, 3, 1],
                [0, 2, 10, 8],
                [9, 11, 3, 1],
                [9, 11, 10, 8],
                [2, 6, 7, 3],
                [2, 6, 14, 10],
                [11, 15, 7, 3],
                [11, 15, 14, 10],
                [0, 4, 5, 1],
                [0, 4, 12, 8],
                [9, 13, 5, 1],
                [9, 13, 12, 8],
                [4, 6, 7, 5],
                [4, 6, 14, 12],
                [13, 15, 14, 12],
                [13, 15, 7, 5],
            ]
        for i in range(nr_treads) :
            if tread_type == TREADTYPE.CLASSIC :
                mm.make_ppd_mesh(treads_coords, 'treads')
            elif tread_type == TREADTYPE.BASIC_STEEL :
                temp = []
                for j in treads_coords :
                    temp.append(copy(j))
                #end for
                for j in range(nr_tread_sections) :
                    mm.make_mesh(temp, basic_steel_faces, 'treads')
                    for k in temp :
                        k += vec(depth + tread_section_spacing, 0, 0)
                    #end for
                #end for
            elif tread_type in [TREADTYPE.BAR_1, TREADTYPE.BAR_2, TREADTYPE.BAR_3] :
                mm.make_mesh(treads_coords, bar_faces, 'treads')
                temp = []
                for j in bars_coords :
                    temp.append(copy(j))
                #end for
                for j in range(nr_tread_sections) :
                    mm.make_ppd_mesh(temp, 'bars')
                    for k in temp :
                        k += vec(tread_section_offset + tread_section_spacing, 0, 0)
                    #end for
                #end for
                for j in bars_coords :
                    j += vec(tread_run, 0, tread_rise)
                #end for
                temp = []
                for j in crosses_coords :
                    temp.append(copy(j))
                #end for
                for j in range(nr_cross_sections) :
                    mm.make_ppd_mesh(temp, 'crosses')
                    for k in temp :
                        k += vec(0, cW + cross, 0)
                    #end for
                #end for
                for j in crosses_coords :
                    j += vec(tread_run, 0, tread_rise)
                #end for
            #end if
            for j in treads_coords :
                j += vec(tread_run, 0, tread_rise)
            #end for
        #end for
    #end make_treads

#begin treads
    make_treads_coords()
    make_treads()
#end treads

def treads_circular(mm, tread_type, stair_arc, outer_radius, tread_height, tread_rise, tread_toe, inner_radius, nr_treads, nr_sections_per_slice = 4) :
    "generates treads for circular stairs."
    start = \
        [
            vec(0, - inner_radius, 0),
            vec(0, - inner_radius, - tread_height),
            vec(0, - outer_radius, 0),
            vec(0, - outer_radius, - tread_height),
        ]
    tread_arc = stair_arc / nr_treads
    for i in range(nr_treads) :
        coords = []
        # Base faces.  Should be able to append more sections:
        bar_2_faces = [[0, 1, 3, 2]]
        t_inner = Matrix.Rotation(- tread_toe / inner_radius + tread_arc * i, 3, 'Z')
        coords.append(t_inner * start[0] + vec(0, 0, tread_rise * i))
        coords.append(t_inner * start[1] + vec(0, 0, tread_rise * i))
        t_outer = Matrix.Rotation(- tread_toe / outer_radius + tread_arc * i, 3, 'Z')
        coords.append(t_outer * start[2] + vec(0, 0, tread_rise * i))
        coords.append(t_outer * start[3] + vec(0, 0, tread_rise * i))
        for j in range(nr_sections_per_slice + 1) :
            k = j * 4 + 4
            bar_2_faces.append([k, k - 4, k - 3, k + 1])
            bar_2_faces.append([k - 2, k - 1, k + 3, k + 2])
            bar_2_faces.append([k + 1, k - 3, k - 1, k + 3])
            bar_2_faces.append([k, k - 4, k - 2, k + 2])
            rot = Matrix.Rotation(tread_arc * j / nr_sections_per_slice + tread_arc * i, 3, 'Z')
            for v in start :
                coords.append(rot * v + vec(0, 0, tread_rise * i))
            #end for
        #end for
        bar_2_faces.append([k, k + 1, k + 3, k + 2])
        mm.make_mesh(coords, bar_2_faces, 'treads')
    #end for
#end treads_circular

#+
# Putting it all together
#-

class Stairs(bpy.types.Operator) :
    "actual generation of stairs."

    bl_idname = "mesh.stairs"
    bl_label = "Add Stairs"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Add stairs"

    stair_type = EnumProperty \
      (
        name = "Type",
        description = "Type of staircase to generate",
        items = STAIRTYPE.all_items()
      )

    rise = FloatProperty \
      (
        name = "Rise",
        description = "Single tread rise",
        min = 0.0,
        max = 1024.0,
        default = 0.20
      )
    run = FloatProperty \
      (
        name = "Run",
        description = "Single tread run",
        min = 0.0,
        max = 1024.0,
        default = 0.30
      )

    #for circular
    rad1 = FloatProperty \
      (
        name = "Inner Radius",
        description = "Inner radius for circular staircase",
        min = 0.0,
        max = 1024.0,
        soft_max = 10.0,
        default = 0.25
      )
    rad2 = FloatProperty \
      (
        name = "Outer Radius",
        description = "Outer radius for circular staircase",
        min = 0.0,
        max = 1024.0,
        soft_min = 0.015625,
        soft_max = 32.0,
        default = 1.0
      )
    rotation = FloatProperty \
      (
        subtype = "ANGLE",
        name = "Rotation",
        description = "How much the stairway rotates",
        min = 0.0,
        max = 256 * circle,
        step = 5 * 100,
        default = 1.25 * circle
      )
    center = BoolProperty \
      (
        name = "Center Pillar",
        description = "Generate a central pillar",
        default = False
      )

    #for treads
    make_treads = BoolProperty \
      (
        name = "Make Treads",
        description = "Enable tread generation",
        default = True
      )
    tread_w = FloatProperty \
      (
        name = "Tread Width",
        description = "Width of each generated tread",
        min = 0.0001,
        max = 1024.0,
        default = 1.2
      )
    tread_h = FloatProperty \
      (
        name = "Tread Height",
        description = "Height of each generated tread",
        min = 0.0001,
        max = 1024.0,
        default = 0.04
      )
    tread_t = FloatProperty \
      (
        name = "Tread Toe",
        description = "Toe (aka \"nosing\") of each generated tread",
        min = 0.0,
        max = 10.0,
        default = 0.03
      )
    tread_o = FloatProperty \
      (
        name = "Tread Overhang",
        description = "How much tread \"overhangs\" the sides",
        min = 0.0,
        max = 1024.0,
        default = 0.025
      )
    tread_n = IntProperty \
      (
        name = "Number of Treads",
        description = "How many treads to generate",
        min = 1,
        max = 1024,
        default = 10
      )
    tread_type = EnumProperty \
      (
        name = "Tread Type",
        description = "Type/style of treads to generate",
        items = TREADTYPE.all_items()
      )
    tread_tk = FloatProperty \
      (
        name = "Thickness",
        description = "Thickness of the treads",
        min = 0.0001, max = 10.0,
        default = 0.02
      )
    tread_sec = IntProperty \
      (
        name = "Sections",
        description = "Number of sections to use for tread",
        min = 1,
        max = 1024,
        default = 5
      )
    tread_sp = IntProperty \
      (
        subtype = "PERCENTAGE",
        name = "Spacing",
        description = "Total spacing between tread sections as a percentage of total tread width",
        min = 0,
        max = 80,
        default = 5
      )
    tread_sn = IntProperty \
      (
        name = "Crosses",
        description = "Number of cross section supports",
        min = 2,
        max = 1024,
        default = 4
      )
    #special circular tread properties:
    tread_slc = IntProperty \
      (
        name = "Slices",
        description = "Number of slices each tread is composed of",
        min = 1,
        max = 1024,
        soft_max = 16,
        default = 4
      )

    #for posts
    make_posts = BoolProperty \
      (
        name = "Make Posts",
        description = "Enable post generation",
        default = True
      )
    post_d = FloatProperty \
      (
        name = "Post Depth",
        description = "Depth of generated posts",
        min = 0.0001,
        max = 10.0,
        default = 0.04
      )
    post_w = FloatProperty \
      (
        name = "Post Width",
        description = "Width of generated posts",
        min = 0.0001,
        max = 10.0,
        default = 0.04
      )
    post_n = IntProperty \
      (
        name = "Number of Posts",
        description = "Number of posts to generated",
        min = 1,
        max = 1024,
        default = 5
      )

    #for railings
    make_railings = BoolProperty \
      (
        name = "Make Railings",
        description = "Generate railings",
        default = True
      )
    rail_w = FloatProperty \
      (
        name = "Railings Width",
        description = "Width of railings to generate",
        min = 0.0001,
        max = 10.0,
        default = 0.12
      )
    rail_t = FloatProperty \
      (
        name = "Railings Thickness",
        description = "Thickness of railings to generate",
        min = 0.0001,
        max = 10.0,
        default = 0.03
      )
    rail_h = FloatProperty \
      (
        name = "Railings Height",
        description = "Height of railings to generate",
        min = 0.0001,
        max = 10.0,
        default = 0.90
      )

    #for retainers
    make_retainers = BoolProperty \
      (
        name = "Make Retainers",
        description = "Generate retainers",
        default = True
      )
    ret_w = FloatProperty \
      (
        name = "Retainer Width",
        description = "Width of generated retainers",
        min = 0.0001,
        max = 10.0,
        default = 0.01
      )
    ret_h = FloatProperty \
      (
        name = "Retainer Height",
        description = "Height of generated retainers",
        min = 0.0001,
        max = 10.0,
        default = 0.01
      )
    ret_n = IntProperty \
      (
        name = "Number of Retainers",
        description = "Number of retainers to generated",
        min = 1,
        max = 1024,
        default = 3
      )

    #for stringer
    make_stringer = BoolProperty \
      (
        name = "Make Stringer",
        description = "Generate stair stringer",
        default = True
      )
    stringer_type = EnumProperty \
      (
        name = "Stringer Type",
        description = "Type/style of stringer to generate",
        items = STRINGERTYPE.all_items()
      )
    string_n = IntProperty \
      (
        name = "Number of Stringers",
        description = "Number of stringers to generate",
        min = 1,
        max = 10,
        default = 1
      )
    string_dis = BoolProperty \
      (
        name = "Distributed",
        description = "Use distributed stringers",
        default = False
      )
    string_w = FloatProperty \
      (
        subtype = "PERCENTAGE",
        name = "Stringer width",
        description = "Width of stringer as a percentage of tread width",
        min = 0.0001,
        max = 100.0,
        default = 15.0
      )
    string_h = FloatProperty \
      (
        name = "Stringer Height",
        description = "Height of the stringer",
        min = 0.0001,
        max = 100.0,
        default = 0.3
      )
    string_tw = FloatProperty \
      (
        subtype = "PERCENTAGE",
        name = "Web Thickness",
        description = "Thickness of the beam's web as a percentage of width",
        min = 0.0001,
        max = 100.0,
        default = 25.0
      )
    string_tf = FloatProperty \
      (
        name = "Flange Thickness",
        description = "Thickness of the flange",
        min = 0.0001,
        max = 100.0,
        default = 0.05
      )
    string_tp = FloatProperty \
      (
        subtype = "PERCENTAGE",
        name = "Flange Taper",
        description = "Flange thickness taper as a percentage",
        min = 0.0,
        max = 100.0,
        default = 0.0
      )
    string_g = BoolProperty \
      (
        name = "Floating",
        description = "Cut bottom of stringer to be a \"floating\" section",
        default = False
      )

    use_original = BoolProperty \
      (
        name = "Use legacy method",
        description = "Use the Blender 2.49 legacy method for stair generation",
        default = True
      )
    rEnable = BoolProperty \
      (
        name = "Right Details",
        description = "Generate right side details (posts/rails/retainers)",
        default = True
      )
    lEnable = BoolProperty \
      (
        name = "Left Details",
        description = "Generate left side details (posts/rails/retainers)",
        default = True
      )

    # Draw the GUI:
    def draw(self, context) :
        layout = self.layout
        box = layout.box()
        box.prop(self, 'stair_type')
        box = layout.box()
        box.prop(self, 'rise')
        if self.stair_type != STAIRTYPE.CIRCULAR.name :
            box.prop(self, 'run')
        else :
            box.prop(self, 'rotation')
            box.prop(self, 'rad1')
            box.prop(self, 'rad2')
            box.prop(self, 'center') # FIXME: not used
        #end if
        if self.stair_type == STAIRTYPE.FREESTANDING.name :
            box.prop(self, 'use_original')
            if not self.use_original :
                box.prop(self, 'rEnable')
                box.prop(self, 'lEnable')
            #end if
        else :
            self.use_original = False
            box.prop(self, 'rEnable')
            box.prop(self, 'lEnable')
        #end if
        # Treads
        box = layout.box()
        box.prop(self, 'make_treads')
        if self.make_treads :
            if not self.use_original and self.stair_type != STAIRTYPE.CIRCULAR.name :
                box.prop(self, 'tread_type')
            else :
                self.tread_type = TREADTYPE.CLASSIC.name
            #end if
            if self.stair_type != STAIRTYPE.CIRCULAR.name :
                box.prop(self, 'tread_w')
            #end if
            box.prop(self, 'tread_h')
            box.prop(self, 'tread_t')
            if self.stair_type not in [STAIRTYPE.HOUSED_OPEN.name, STAIRTYPE.CIRCULAR.name] :
                box.prop(self, 'tread_o')
            else :
                self.tread_o = 0.0
            #end if
            box.prop(self, 'tread_n')
            if self.tread_type != TREADTYPE.CLASSIC.name :
                box.prop(self, 'tread_tk')
                box.prop(self, 'tread_sec')
                if self.tread_sec > 1 and self.tread_type not in [TREADTYPE.BAR_1.name, TREADTYPE.BAR_2.name] :
                    box.prop(self, 'tread_sp')
                #end if
                if self.tread_type in [TREADTYPE.BAR_1.name, TREADTYPE.BAR_2.name, TREADTYPE.BAR_3.name] :
                    box.prop(self, 'tread_sn')
                #end if
            elif self.stair_type == STAIRTYPE.CIRCULAR.name :
                box.prop(self, "tread_slc")
            #end if
        #end if
        # Posts
        box = layout.box()
        box.prop(self, 'make_posts')
        if self.make_posts :
            box.prop(self, 'post_d')
            box.prop(self, 'post_w')
            box.prop(self, 'post_n')
        #end if
        # Railings
        box = layout.box()
        box.prop(self, 'make_railings')
        if self.make_railings :
            box.prop(self, 'rail_w')
            box.prop(self, 'rail_t')
            box.prop(self, 'rail_h')
        #end if
        # Retainers
        box = layout.box()
        box.prop(self, 'make_retainers')
        if self.make_retainers :
            box.prop(self, 'ret_w')
            box.prop(self, 'ret_h')
            box.prop(self, 'ret_n')
        #end if
        # Stringers
        box = layout.box()
        if self.stair_type != STAIRTYPE.HOUSED_OPEN.name :
            box.prop(self, 'make_stringer')
        else :
            self.make_stringer = True
        #end if
        if self.make_stringer :
            if not self.use_original :
                box.prop(self, 'stringer_type')
            else :
                self.stringer_type = STRINGERTYPE.CLASSIC.name
            #end if
            box.prop(self, 'string_w')
            if self.stair_type == STAIRTYPE.FREESTANDING.name :
                if self.stringer_type == STRINGERTYPE.CLASSIC.name and not self.use_original :
                    box.prop(self, 'string_n')
                    box.prop(self, 'string_dis')
                elif self.stringer_type in [STRINGERTYPE.I_BEAM.name, STRINGERTYPE.C_BEAM.name] :
                    box.prop(self, 'string_n')
                    box.prop(self, 'string_dis')
                    box.prop(self, 'string_h')
                    box.prop(self, 'string_tw')
                    box.prop(self, 'string_tf')
                    box.prop(self, 'string_tp')
                    box.prop(self, 'string_g')
                #end if
            elif self.stair_type == STAIRTYPE.HOUSED_OPEN.name :
                if self.stringer_type in [STRINGERTYPE.I_BEAM.name, STRINGERTYPE.C_BEAM.name] :
                    box.prop(self, 'string_tw')
                    box.prop(self, 'string_tf')
                #end if
            #end if
        #end if
        # Tread support:
##        if self.make_stringer and self.stringer_type in [STRINGERTYPE.I_BEAM.name, STRINGERTYPE.C_BEAM.name] :
    #end draw

    def execute(self, context) :
        mm = MeshMaker(self.rise, self.run, self.tread_n)
        # convert strings to enums
        stair_type = STAIRTYPE[self.stair_type]
        stringer_type = STRINGERTYPE[self.stringer_type]
        tread_type = TREADTYPE[self.tread_type]
        if self.make_treads :
            if stair_type != STAIRTYPE.CIRCULAR :
                treads \
                  (
                    mm = mm,
                    stair_type = stair_type,
                    tread_type = tread_type,
                    stair_run = self.run,
                    tread_width = self.tread_w,
                    tread_height = self.tread_h,
                    tread_run = self.run,
                    tread_rise = self.rise,
                    tread_toe = self.tread_t,
                    tread_side_overhang = self.tread_o,
                    nr_treads = self.tread_n,
                    tread_metal_thickness = self.tread_tk,
                    nr_tread_sections = self.tread_sec,
                    spacing = self.tread_sp,
                    nr_cross_sections = self.tread_sn
                  )
            else :
                treads_circular \
                  (
                    mm = mm,
                    tread_type = tread_type,
                    stair_arc = self.rotation,
                    outer_radius = self.rad2,
                    tread_height = self.tread_h,
                    tread_rise = self.rise,
                    tread_toe = self.tread_t,
                    inner_radius = self.rad1,
                    nr_treads = self.tread_n,
                    nr_sections_per_slice = self.tread_slc
                  )
            #end if
        #end if
        if self.make_posts and (self.rEnable or self.lEnable) :
            posts \
              (
                mm = mm,
                rise = self.rise,
                stair_run = self.run,
                post_depth = self.post_d,
                post_width = self.post_w,
                tread_width = self.tread_w,
                nr_posts = self.post_n,
                rail_height = self.rail_h,
                rail_thickness = self.rail_t,
                rEnable = self.rEnable,
                lEnable = self.lEnable
              )
        #end if
        if self.make_railings and (self.rEnable or self.lEnable) :
            railings \
              (
                mm = mm,
                rail_width = self.rail_w,
                rail_thickness = self.rail_t,
                rail_height = self.rail_h,
                tread_toe = self.tread_t,
                post_width = self.post_w,
                post_depth = self.post_d,
                tread_width = self.tread_w,
                rEnable = self.rEnable,
                lEnable = self.lEnable
              )
        #end if
        if self.make_retainers and (self.rEnable or self.lEnable) :
            retainers \
              (
                mm = mm,
                retainer_width = self.ret_w,
                retainer_height = self.ret_h,
                post_width = self.post_w,
                tread_width = self.tread_w,
                rail_height = self.rail_h,
                nr_retainers = self.ret_n,
                rEnable = self.rEnable,
                lEnable = self.lEnable
              )
        #end if
        if self.make_stringer :
            if stair_type == STAIRTYPE.FREESTANDING and self.use_original :
                stringer \
                  (
                    mm = mm,
                    stair_type = stair_type,
                    stringer_type = stringer_type,
                    stair_rise = self.rise,
                    stair_run = self.run,
                    w = self.string_w,
                    stringer_height = self.string_h,
                    nr_treads = self.tread_n,
                    tread_height = self.tread_h,
                    tread_width = self.tread_w,
                    tread_toe = self.tread_t,
                    tread_overhang = self.tread_o,
                    tw = self.string_tw,
                    stringer_flange_thickness = self.string_tf,
                    tp = self.string_tp,
                    stringer_intersects_ground = not self.string_g
                  )
            elif stair_type == STAIRTYPE.BOX :
                stringer \
                  (
                    mm = mm,
                    stair_type = stair_type,
                    stringer_type = stringer_type,
                    stair_rise = self.rise,
                    stair_run = self.run,
                    w = 100,
                    stringer_height = self.string_h,
                    nr_treads = self.tread_n,
                    tread_height = self.tread_h,
                    tread_width = self.tread_w,
                    tread_toe = self.tread_t,
                    tread_overhang = self.tread_o,
                    tw = self.string_tw,
                    stringer_flange_thickness = self.string_tf,
                    tp = self.string_tp,
                    stringer_intersects_ground = not self.string_g,
                    nr_stringers = 1,
                    distributed_stringers = False,
                    notMulti = False
                  )
            elif stair_type == STAIRTYPE.CIRCULAR :
                stringer \
                  (
                    mm = mm,
                    stair_type = stair_type,
                    stringer_type = stringer_type,
                    stair_rise = self.rise,
                    stair_run = self.rotation,
                    w = self.string_w,
                    stringer_height = self.string_h,
                    nr_treads = self.tread_n,
                    tread_height = self.tread_h,
                    tread_width = self.rad2 - self.rad1,
                    tread_toe = self.tread_t,
                    tread_overhang = self.rad1,
                    tw = self.string_tw,
                    stringer_flange_thickness = self.string_tf,
                    tp = self.string_tp,
                    stringer_intersects_ground = not self.string_g,
                    nr_stringers = self.string_n,
                    distributed_stringers = self.string_dis,
                    notMulti = self.use_original,
                    sections_per_slice = self.tread_slc
                  )
            else :
                stringer \
                  (
                    mm = mm,
                    stair_type = stair_type,
                    stringer_type = stringer_type,
                    stair_rise = self.rise,
                    stair_run = self.run,
                    w = self.string_w,
                    stringer_height = self.string_h,
                    nr_treads = self.tread_n,
                    tread_height = self.tread_h,
                    tread_width = self.tread_w,
                    tread_toe = self.tread_t,
                    tread_overhang = self.tread_o,
                    tw = self.string_tw,
                    stringer_flange_thickness = self.string_tf,
                    tp = self.string_tp,
                    stringer_intersects_ground = not self.string_g,
                    nr_stringers = self.string_n,
                    distributed_stringers = self.string_dis,
                    notMulti = self.use_original
                  )
            #end if
        #end if
        if len(mm.made_objects) != 0 :
            # cannot seem to create an empty with bpy.data.objects.new()
            bpy.ops.object.empty_add \
              (
                location = bpy.context.scene.cursor_location,
                layers = bpy.context.space_data.layers
              )
            root = bpy.context.active_object
            assert root.type == "EMPTY"
            root.name = "staircase root"
            for obj in mm.made_objects :
                obj.parent = root
                obj.parent_type = "OBJECT"
                obj.select = True
            #end if
            bpy.ops.group.create(name = "staircase components")
            for obj in mm.made_objects :
                obj.select = False
            #end if
        #end if
        return {'FINISHED'}
    #end execute

#end Stairs
