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
#   - global vs. local needs cleaned up.
#   - Join separate stringer objects and then clean up the mesh.
#   - Put all objects into a group.
#   - Generate left/right posts/railings/retainers separately with
#       option to disable just the left/right.
#   - Add wall railing type as an option for left/right
#   - Add different rail styles (profiles).  Select with enum.
#   - Would like to add additional staircase types.
#       - Spiral staircase
#       - "L" staircase
#       - "T" staircase
#
# Last Modified By: Paul "brikbot" Marshall
# Last Modification: January 29, 2011
#
# ##### BEGIN GPL LICENSE BLOCK #####
#
#  Stairbuilder is for quick stair generation.
#  Copyright (C) 2010  Nick van Adium
#  Copyright (C) 2011  Paul Marshall
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
from copy import copy
import bpy
from bpy.props import (BoolProperty,
                       EnumProperty,
                       IntProperty,
                       FloatProperty)
from bpy_extras import object_utils
from mathutils import \
    Matrix, \
    Vector
from mathutils.geometry import \
    intersect_line_plane, \
    intersect_line_line

global G
global typ
global typ_s
global typ_t
global rise
global run

# from general.py

class General:
    def __init__(self,rise,run,N):
        self.stop=float(N)*Vector([run,0,rise])
        self.slope=rise/run
        self.angle=math.atan(self.slope)
        #identical quads for all objects except stringer
        self.faces=[[0,1,3,2],[0,1,5,4],[0,2,6,4],[4,5,7,6],[2,3,7,6],[1,3,7,5]]

    def Make_mesh(self, verts, faces, name):
        # Create new mesh
        mesh = bpy.data.meshes.new(name)

        # Make a mesh from a list of verts/edges/faces.
        mesh.from_pydata(verts, [], faces)

        # Set mesh to use auto smoothing:
        mesh.use_auto_smooth = True

        # Update mesh geometry after adding stuff.
        mesh.update()

        return object_utils.object_data_add(bpy.context, mesh, operator=None)

# from post.py

class Posts:
    def __init__(self,G,rise,run,d,w,wT,nP,hR,tR, rEnable, lEnable):
        self.G = G #General
        self.rise = rise #Stair rise
        self.run = run #Stair run
        self.x1=Vector([0,0,hR-tR]) #rail start
        self.x2=G.stop+Vector([0,0,hR-tR]) #rail stop
        self.d=d #post depth
        self.w=w #post width
        self.wT=wT #tread width
        self.nP=nP #number of posts
        self.sp=Vector([(self.x2[0]-self.x1[0])/float(nP+1),0,0]) #spacing between posts
        self.rEnable = rEnable
        self.lEnable = lEnable
        self.Create()

    def Intersect(self,i,d):
        """find intersection point, x, for rail and post"""
        x3=self.x1+i*self.sp+Vector([d,d,d])
        x4=x3+Vector([0,0,self.x2[-1]])
        a=self.x2-self.x1
        b=x4-x3
        c=x3-self.x1
        cr_ab=a.cross(b)
        mag_cr_ab=(cr_ab * cr_ab)
        return self.x1+a*((c.cross(b).dot(cr_ab))/mag_cr_ab)

    def Create(self):
        for i in range(0,self.nP+2,1):
            coords = []
            #intersections with rail
            coords.append(self.Intersect(i,0.0))
            coords.append(self.Intersect(i,self.d))
            #intersections with tread
            coords.append(Vector([self.x1[0]+i*self.sp[0],0,
                                  int(coords[0][0]/self.run)*self.rise]))
            coords.append(coords[2]+Vector([self.d,0,0]))
            #inner face
            for j in range(4):
                coords.append(coords[j]+Vector([0,self.w,0]))
            if self.rEnable:
                self.G.Make_mesh(coords, self.G.faces, 'posts')
            if self.lEnable:
                #make post on other side of steps as well
                for j in coords:
                    j += Vector([0,self.wT-self.w,0])
                self.G.Make_mesh(coords, self.G.faces, 'posts')

# from rail.py

class Rails:
    def __init__(self,G,w,t,h,tT,wP,dP,wT, rEnable, lEnable):
        self.G = G #General
        self.w=w #rail width
        self.t=t #rail thickness
        self.h=h #rail height
        self.start=Vector([0,0,self.h-self.t]) #rail start
        self.stop=G.stop+Vector([0,0,self.h-self.t]) #rail stop
        self.tT=tT #tread toe
        self.wP=wP #post width
        self.dP=dP #post depth
        self.wT=wT #tread width
        self.rEnable = rEnable
        self.lEnable = lEnable
        self.Create()

    def Create(self):
        #determine offset to include railing toe
        offset=Vector([self.tT,0,self.tT*math.tan(self.G.angle)])
        coords = []
        coords.append(self.start-offset)
        coords.append(self.stop+offset+Vector([self.dP,0,
                                               self.dP*math.tan(self.G.angle)]))
        coords.append(self.start-offset+Vector([0,self.w,0]))
        coords.append(self.stop+offset+Vector([self.dP,self.w,
                                               self.dP*math.tan(self.G.angle)]))
        for j in range(4):
            coords.append(coords[j]+Vector([0,0,self.t]))
        #centre over posts
        for j in coords:
            j += Vector([0,0.5*(-self.w+self.wP),0])
        if self.rEnable:
            self.G.Make_mesh(coords, self.G.faces, 'rails')
        if self.lEnable:
            #make rail on other side
            for j in coords:
                j += Vector([0,self.wT-self.wP,0])
            self.G.Make_mesh(coords, self.G.faces, 'rails')

# from retainer.py

class Retainers:
    def __init__(self,G,w,h,wP,wT,hR,n, rEnable, lEnable):
        self.G = G #General
        self.w=w #retainer width
        self.h=h #retainer height
        self.wP=wP #post width
        self.wT=wT #tread width
        self.nR=n #number of retainers
        self.sp=hR/float(n+1) #retainer spacing
        self.rEnable = rEnable
        self.lEnable = lEnable
        self.Create()

    def Create(self):
        for i in range(self.nR):
            coords = []
            offset=(i+1)*Vector([0,0,self.sp])
            coords.append(offset)
            coords.append(self.G.stop + offset)
            coords.append(offset + Vector([0,self.w,0]))
            coords.append(self.G.stop + offset + Vector([0,self.w,0]))
            for j in range(4):
                coords.append(coords[j] + Vector([0,0,self.h]))
            #centre in posts
            for j in coords:
                j += Vector([0,0.5*(self.wP-self.w),0])
            if self.rEnable:
                self.G.Make_mesh(coords, self.G.faces, 'retainers')
            if self.lEnable:
                #make retainer on other side
                for j in coords:
                    j += Vector([0,self.wT-self.wP,0])
                self.G.Make_mesh(coords,self.G.faces, 'retainers')

# from stringer.py

class Stringer:
    def  __init__(self,G,typ,typ_s,rise,run,w,h,nT,hT,wT,tT,tO,tw,tf,tp,g,
                  nS=1,dis=False,notMulti=True,deg=4):
        self.G = G #General
        self.typ = typ # Stair type
        self.typ_s = typ_s # Stringer type
        self.rise = rise #Stair rise
        self.run = run #Stair run. Degrees if self.typ == "id4"
        if notMulti:
            self.w = w / 100 #stringer width
        else:
            self.w = (wT * (w / 100)) / nS
        self.h = h #stringer height
        self.nT = nT #number of treads
        self.hT = hT #tread height
        self.wT = wT #tread width
        self.tT = tT #tread toe
        self.tO = tO #Tread overhang. Inner radius if self.typ == "id4"
        self.tw = self.w * (tw / 100) #stringer web thickness
        self.tf = tf #stringer flange thickness
        self.tp = 1 - (tp / 100) #stringer flange taper
        self.g = g #does stringer intersect the ground?
        self.nS = nS #number of stringers
        self.dis = dis #Use distributed stringers
        self.deg = deg #number of sections per "slice". Only applys if self.typ == "id4"
        # Default stringer object (classic / sId1):
        self.faces1=[[0,1,3,2],[1,5,3],[3,5,4],[6,7,9,8],[7,11,9],[9,11,10],
                     [0,2,8,6],[0,1,7,6],[1,5,11,7],[2,3,9,8],[3,4,10,9],[4,5,11,10]]
        # Box stair type stringer:
        self.faces2=[[0,1,7,6],[1,3,9,7],[3,4,10,9],[4,10,11,5],[5,11,8,2],
                     [2,8,6,0],[0,1,2],[1,2,5,3],[3,4,5],[6,7,8],[7,8,11,9],[9,10,11]]
        # I-beam stringer (id2 / sId2 / Taper < 100%):
        self.faces3a=[[0,1,17,16],[1,2,18,17],[2,3,19,18],[3,4,20,19],[4,5,21,20],[5,6,22,21],
                      [6,7,23,22],[7,8,24,23],[8,9,25,24],[9,10,26,25],[10,11,27,26],
                      [11,12,28,27],[12,13,29,28],[13,14,30,29],[14,15,31,30],[15,0,16,31],
                      [0,1,2,15],[2,11,14,15],[11,12,13,14],[2,3,10,11],[3,4,5,6],[3,6,7,10],
                      [7,8,9,10],[16,17,18,31],[18,27,30,31],[27,28,29,30],[18,19,26,27],
                      [19,20,21,22],[19,22,23,26],[23,24,25,26]]
        # I-beam stringer (id2 / sId2 / Taper = 100%):
        self.faces3b=[[0,1,9,8],[1,2,10,9],[2,3,11,10],[3,4,12,11],[4,5,13,12],[5,6,14,13],
                      [6,7,15,14],[7,0,8,15],[0,1,6,7],[1,2,5,6],[2,3,4,5],[8,9,14,15],
                      [9,10,13,14],[10,11,12,13]]
        # I-beam stringer (id3 / sId2 / Taper < 100%):
        self.faces3c=[[0,1,2,7],[2,3,6,7],[3,4,5,6],[1,2,23,16],[2,3,22,23],
                      [3,4,21,22],[16,17,18,23],[18,19,22,23],[19,20,21,22],
                      [17,8,15,18],[18,15,14,19],[19,14,13,20],[8,9,10,15],
                      [10,11,14,15],[11,12,13,14],[9,10,53,52],[10,11,54,53],
                      [11,12,55,54],[52,53,61,60],[53,54,62,61],[54,55,63,62],
                      [60,61,34,33],[61,62,35,34],[62,63,36,35],[32,33,34,39],
                      [34,35,38,39],[35,36,37,38],[41,32,39,42],[42,39,38,43],
                      [43,38,37,44],[40,41,42,47],[42,43,46,47],[43,44,45,46],
                      [25,26,47,40],[26,27,46,47],[27,28,45,46],[24,25,26,31],
                      [26,27,30,31],[27,28,29,30],[24,31,57,56],[31,30,58,57],
                      [30,29,59,58],[48,49,57,56],[49,50,58,57],[50,51,59,58],
                      [0,7,49,48],[7,6,50,49],[6,5,51,50],[0,1,16,48],[16,40,56,48],
                      [24,25,40,56],[16,17,41,40],[8,9,52,17],[17,52,60,41],
                      [32,33,60,41],[12,13,20,55],[20,44,63,55],[37,44,63,36],
                      [20,21,45,44],[28,29,51,21],[21,51,59,45],[28,45,59,29],
                      [4,5,51,21]]
        # C-beam stringer (id3 / sId3 / Taper < 100%):
        self.faces4c=[[0,1,2,7],[2,3,6,7],[3,4,5,6],[1,2,23,16],[2,3,22,23],[3,4,21,22],
                      [16,17,18,23],[18,19,22,23],[19,20,21,22],[17,8,15,18],[18,15,14,19],
                      [19,14,13,20],[8,9,10,15],[10,11,14,15],[11,12,13,14],[0,24,25,7],
                      [7,25,26,6],[6,26,27,5],[9,31,30,10],[10,30,29,11],[11,29,28,12],
                      [24,25,30,31],[25,26,29,30],[26,27,28,29],[0,1,16,24],[16,24,31,17],
                      [8,9,31,17],[4,5,27,21],[20,21,27,28],[12,13,20,28]]
        self.Create()


    def Create(self):
        if self.typ == "id1":
            if self.typ_s == "sId1":
                if self.dis or self.nS == 1:
                    offset = (self.wT / (self.nS + 1)) - (self.w / 2)
                else:
                    offset = 0
                for i in range(self.nS):
                    for j in range(self.nT):
                        coords = []
                        coords.append(Vector([0, offset, -self.rise]))
                        coords.append(Vector([self.run, offset, -self.rise]))
                        coords.append(Vector([0, offset, -self.hT]))
                        coords.append(Vector([self.run, offset, -self.hT]))
                        coords.append(Vector([self.run, offset, 0]))
                        coords.append(Vector([self.run * 2, offset, 0]))
                        for k in range(6):
                            coords.append(coords[k]+Vector([0, self.w, 0]))
                        for k in coords:
                            k += j*Vector([self.run, 0, self.rise])
                        self.G.Make_mesh(coords,self.faces1,'stringer')
                    if self.dis or self.nS == 1:
                        offset += self.wT / (self.nS + 1)
                    else:
                        offset += (self.wT - self.w) / (self.nS - 1)
            elif self.typ_s == "sId2":
                self.I_beam()
        elif self.typ == "id2":
            if self.typ_s == "sId1":
                coords = []
                coords.append(Vector([-self.tT, -self.w, -self.rise]))
                coords.append(Vector([self.hT / self.G.slope, -self.w, -self.rise]))
                coords.append(Vector([-self.tT, -self.w, 0]))
                coords.append(Vector([self.nT * self.run, -self.w,
                                      ((self.nT - 1) * self.rise) - self.hT]))
                coords.append(Vector([self.nT * self.run, -self.w, self.nT * self.rise]))
                coords.append(Vector([(self.nT * self.run) - self.tT, -self.w,
                                      self.nT * self.rise]))
                for i in range(6):
                    coords.append(coords[i] + Vector([0, self.w, 0]))
                self.G.Make_mesh(coords, self.faces2, 'stringer')
                for i in coords:
                    i += Vector([0, self.w + self.wT, 0])
                self.G.Make_mesh(coords, self.faces2, 'stringer')
            elif self.typ_s == "sId2":
                self.housed_I_beam()
            elif self.typ_s == "sId3":
                self.housed_C_beam()
        elif self.typ == "id3":
            h = (self.rise - self.hT) - self.rise #height of top section
            for i in range(self.nT):
                coords = []
                coords.append(Vector([i * self.run,0,-self.rise]))
                coords.append(Vector([(i + 1) * self.run,0,-self.rise]))
                coords.append(Vector([i * self.run,0,h + (i * self.rise)]))
                coords.append(Vector([(i + 1) * self.run,0,h + (i * self.rise)]))
                for j in range(4):
                    coords.append(coords[j] + Vector([0,self.wT,0]))
                self.G.Make_mesh(coords, self.G.faces, 'stringer')
        elif self.typ == "id4":
            offset = (self.wT / (self.nS + 1)) - (self.w / 2)
            for s in range(self.nS):
                base = self.tO + (offset * (s + 1))
                start = [Vector([0, -base, -self.hT]),
                         Vector([0, -base, -self.hT - self.rise]),
                         Vector([0, -base - self.w, -self.hT]),
                         Vector([0, -base - self.w, -self.hT - self.rise])]
                self.d = math.radians(self.run) / self.nT
                for i in range(self.nT):
                    coords = []
                    # Base faces.  Should be able to append more sections:
                    tId4_faces = [[0, 1, 3, 2]]
                    t_inner = Matrix.Rotation(self.d * i, 3, 'Z')
                    coords.append((t_inner * start[0]) + Vector([0, 0, self.rise * i]))
                    coords.append((t_inner * start[1]) + Vector([0, 0, self.rise * i]))
                    t_outer = Matrix.Rotation(self.d * i, 3, 'Z')
                    coords.append((t_outer * start[2]) + Vector([0, 0, self.rise * i]))
                    coords.append((t_outer * start[3]) + Vector([0, 0, self.rise * i]))
                    k = 0
                    for j in range(self.deg):
                        k = (j * 4) + 4
                        tId4_faces.append([k, k - 4, k - 3, k + 1])
                        tId4_faces.append([k - 2, k - 1, k + 3, k + 2])
                        tId4_faces.append([k + 1, k - 3, k - 1, k + 3])
                        tId4_faces.append([k, k - 4, k - 2, k + 2])
                        rot = Matrix.Rotation(((self.d * (j + 1)) / self.deg) + (self.d * i), 3, 'Z')
                        for v in start:
                            coords.append((rot * v) + Vector([0, 0, self.rise * i]))
                    for j in range(self.deg):
                        k = ((j + self.deg) * 4) + 4
                        tId4_faces.append([k, k - 4, k - 3, k + 1])
                        tId4_faces.append([k - 2, k - 1, k + 3, k + 2])
                        tId4_faces.append([k + 1, k - 3, k - 1, k + 3])
                        tId4_faces.append([k, k - 4, k - 2, k + 2])
                        rot = Matrix.Rotation(((self.d * ((j + self.deg) + 1)) / self.deg) + (self.d * i), 3, 'Z')
                        for v in range(4):
                            if v in [1, 3]:
                                incline = (self.rise * i) + (self.rise / self.deg) * (j + 1)
                                coords.append((rot * start[v]) + Vector([0, 0, incline]))
                            else:
                                coords.append((rot * start[v]) + Vector([0, 0, self.rise * i]))
                    self.G.Make_mesh(coords, tId4_faces, 'treads')

        return {'FINISHED'}


    def I_beam(self):
        mid = self.w / 2
        web = self.tw / 2
        # Bottom of the stringer:
        baseZ = -self.rise - self.hT - self.h
        # Top of the strigner:
        topZ = -self.rise - self.hT
        # Vertical taper amount:
        taper = self.tf * self.tp

        if self.dis or self.nS == 1:
            offset = (self.wT / (self.nS + 1)) - mid
        else:
            offset = 0

        # taper < 100%:
        if self.tp > 0:
            for i in range(self.nS):
                coords = []
                coords.append(Vector([0, offset,                baseZ]))
                coords.append(Vector([0, offset,                baseZ + taper]))
                coords.append(Vector([0, offset + (mid - web),  baseZ + self.tf]))
                coords.append(Vector([0, offset + (mid - web),  topZ - self.tf]))
                coords.append(Vector([0, offset,                topZ - taper]))
                coords.append(Vector([0, offset,                topZ]))
                coords.append(Vector([0, offset + (mid - web),  topZ]))
                coords.append(Vector([0, offset + (mid + web),  topZ]))
                coords.append(Vector([0, offset + self.w,       topZ]))
                coords.append(Vector([0, offset + self.w,       topZ - taper]))
                coords.append(Vector([0, offset + (mid + web),  topZ - self.tf]))
                coords.append(Vector([0, offset + (mid + web),  baseZ + self.tf]))
                coords.append(Vector([0, offset + self.w,       baseZ + taper]))
                coords.append(Vector([0, offset + self.w,       baseZ]))
                coords.append(Vector([0, offset + (mid + web),  baseZ]))
                coords.append(Vector([0, offset + (mid - web),  baseZ]))
                for j in range(16):
                    coords.append(coords[j]+Vector([self.run * self.nT, 0, self.rise * self.nT]))
                # If the bottom meets the ground:
                #   Bottom be flat with the xy plane, but shifted down.
                #   Either project onto the plane along a vector (hard) or use the built in
                #       interest found in mathutils.geometry (easy).  Using intersect:
                if self.g:
                    for j in range(16):
                        coords[j] = intersect_line_plane(coords[j], coords[j + 16],
                                                         Vector([0, 0, topZ]),
                                                         Vector([0, 0, 1]))
                self.G.Make_mesh(coords, self.faces3a, 'stringer')

                if self.dis or self.nS == 1:
                    offset += self.wT / (self.nS + 1)
                else:
                    offset += (self.wT - self.w) / (self.nS - 1)
        # taper = 100%:
        else:
            for i in range(self.nS):
                coords = []
                coords.append(Vector([0, offset,                baseZ]))
                coords.append(Vector([0, offset + (mid - web),  baseZ + self.tf]))
                coords.append(Vector([0, offset + (mid - web),  topZ - self.tf]))
                coords.append(Vector([0, offset,                topZ]))
                coords.append(Vector([0, offset + self.w,       topZ]))
                coords.append(Vector([0, offset + (mid + web),  topZ - self.tf]))
                coords.append(Vector([0, offset + (mid + web),  baseZ + self.tf]))
                coords.append(Vector([0, offset + self.w,       baseZ]))
                for j in range(8):
                    coords.append(coords[j]+Vector([self.run * self.nT, 0, self.rise * self.nT]))
                self.G.Make_mesh(coords, self.faces3b, 'stringer')
                offset += self.wT / (self.nS + 1)

        return {'FINISHED'}


    def housed_I_beam(self):
        webOrth = Vector([self.rise, 0, -self.run]).normalized()
        webHeight = Vector([self.run + self.tT, 0, -self.hT]).project(webOrth).length
        vDelta_1 = self.tf * math.tan(self.G.angle)
        vDelta_2 = (self.rise * (self.nT - 1)) - (webHeight + self.tf)
        flange_y = (self.w - self.tw) / 2
        front = -self.tT - self.tf
        outer = -self.tO - self.tw - flange_y

        coords = []
        if self.tp > 0:
            # Upper-Outer flange:
            coords.append(Vector([front, outer, -self.rise]))
            coords.append(Vector([-self.tT, outer, -self.rise]))
            coords.append(Vector([-self.tT, outer, 0]))
            coords.append(Vector([(self.run * (self.nT - 1)) - self.tT, outer,
                                  self.rise * (self.nT - 1)]))
            coords.append(Vector([self.run * self.nT, outer,
                                  self.rise * (self.nT - 1)]))
            coords.append(Vector([self.run * self.nT, outer,
                                  (self.rise * (self.nT - 1)) + self.tf]))
            coords.append(Vector([(self.run * (self.nT - 1)) - self.tT, outer,
                                  (self.rise * (self.nT - 1)) + self.tf]))
            coords.append(Vector([front, outer, self.tf - vDelta_1]))
            # Lower-Outer flange:
            coords.append(coords[0] + Vector([self.tf + webHeight, 0, 0]))
            coords.append(coords[1] + Vector([self.tf + webHeight, 0, 0]))
            coords.append(intersect_line_line(coords[9],
                                              coords[9] - Vector([0, 0, 1]),
                                              Vector([self.run, 0, -self.hT - self.tf]),
                                              Vector([self.run * 2, 0, self.rise - self.hT - self.tf]))[0])
            coords.append(Vector([(self.run * self.nT) - ((webHeight - self.hT) / math.tan(self.G.angle)),
                                  outer, vDelta_2]))
            coords.append(coords[4] - Vector([0, 0, self.tf + webHeight]))
            coords.append(coords[5] - Vector([0, 0, self.tf + webHeight]))
            coords.append(coords[11] + Vector([0, 0, self.tf]))
            coords.append(intersect_line_line(coords[8],
                                              coords[8] - Vector([0, 0, 1]),
                                              Vector([self.run, 0, -self.hT]),
                                              Vector([self.run * 2, 0, self.rise - self.hT]))[0])
            # Outer web:
            coords.append(coords[1] + Vector([0, flange_y, 0]))
            coords.append(coords[8] + Vector([0, flange_y, 0]))
            coords.append(coords[15] + Vector([0, flange_y, 0]))
            coords.append(coords[14] + Vector([0, flange_y, 0]))
            coords.append(coords[13] + Vector([0, flange_y, 0]))
            coords.append(coords[4] + Vector([0, flange_y, 0]))
            coords.append(coords[3] + Vector([0, flange_y, 0]))
            coords.append(coords[2] + Vector([0, flange_y, 0]))
            # Upper-Inner flange and lower-inner flange:
            for i in range(16):
                coords.append(coords[i] + Vector([0, self.w, 0]))
            # Inner web:
            for i in range(8):
                coords.append(coords[i + 16] + Vector([0, self.tw, 0]))
            # Mid nodes to so faces will be quads:
            for i in [0,7,6,5,9,10,11,12]:
                coords.append(coords[i] + Vector([0, flange_y, 0]))
            for i in range(8):
                coords.append(coords[i + 48] + Vector([0, self.tw, 0]))

            self.G.Make_mesh(coords, self.faces3c, 'stringer')

            for i in coords:
                i += Vector([0, self.wT + self.tw, 0])

            self.G.Make_mesh(coords, self.faces3c, 'stringer')

        # @TODO Taper = 100%

        return {'FINISHED'}


    def C_Beam(self):
        mid = self.w / 2
        web = self.tw / 2
        # Bottom of the stringer:
        baseZ = -self.rise - self.hT - self.h
        # Top of the strigner:
        topZ = -self.rise - self.hT
        # Vertical taper amount:
        taper = self.tf * self.tp

        if self.dis or self.nS == 1:
            offset = (self.wT / (self.nS + 1)) - mid
        else:
            offset = 0

        # taper < 100%:
        if self.tp > 0:
            for i in range(self.nS):
                coords = []
                coords.append(Vector([0, offset,                baseZ]))
                coords.append(Vector([0, offset,                baseZ + taper]))
                coords.append(Vector([0, offset + (mid - web),  baseZ + self.tf]))
                coords.append(Vector([0, offset + (mid - web),  topZ - self.tf]))
                coords.append(Vector([0, offset,                topZ - taper]))
                coords.append(Vector([0, offset,                topZ]))
                coords.append(Vector([0, offset + (mid - web),  topZ]))
                coords.append(Vector([0, offset + (mid + web),  topZ]))
                coords.append(Vector([0, offset + self.w,       topZ]))
                coords.append(Vector([0, offset + self.w,       topZ - taper]))
                coords.append(Vector([0, offset + (mid + web),  topZ - self.tf]))
                coords.append(Vector([0, offset + (mid + web),  baseZ + self.tf]))
                coords.append(Vector([0, offset + self.w,       baseZ + taper]))
                coords.append(Vector([0, offset + self.w,       baseZ]))
                coords.append(Vector([0, offset + (mid + web),  baseZ]))
                coords.append(Vector([0, offset + (mid - web),  baseZ]))
                for j in range(16):
                    coords.append(coords[j]+Vector([self.run * self.nT, 0, self.rise * self.nT]))
                # If the bottom meets the ground:
                #   Bottom be flat with the xy plane, but shifted down.
                #   Either project onto the plane along a vector (hard) or use the built in
                #       interest found in mathutils.geometry (easy).  Using intersect:
                if self.g:
                    for j in range(16):
                        coords[j] = intersect_line_plane(coords[j], coords[j + 16],
                                                         Vector([0, 0, topZ]),
                                                         Vector([0, 0, 1]))
                self.G.Make_mesh(coords, self.faces3a, 'stringer')

                if self.dis or self.nS == 1:
                    offset += self.wT / (self.nS + 1)
                else:
                    offset += (self.wT - self.w) / (self.nS - 1)
        # taper = 100%:
        else:
            for i in range(self.nS):
                coords = []
                coords.append(Vector([0, offset,                baseZ]))
                coords.append(Vector([0, offset + (mid - web),  baseZ + self.tf]))
                coords.append(Vector([0, offset + (mid - web),  topZ - self.tf]))
                coords.append(Vector([0, offset,                topZ]))
                coords.append(Vector([0, offset + self.w,       topZ]))
                coords.append(Vector([0, offset + (mid + web),  topZ - self.tf]))
                coords.append(Vector([0, offset + (mid + web),  baseZ + self.tf]))
                coords.append(Vector([0, offset + self.w,       baseZ]))
                for j in range(8):
                    coords.append(coords[j]+Vector([self.run * self.nT, 0, self.rise * self.nT]))
                self.G.Make_mesh(coords, self.faces3b, 'stringer')
                offset += self.wT / (self.nS + 1)

        return {'FINISHED'}


    def housed_C_beam(self):
        webOrth = Vector([self.rise, 0, -self.run]).normalized()
        webHeight = Vector([self.run + self.tT, 0, -self.hT]).project(webOrth).length
        vDelta_1 = self.tf * math.tan(self.G.angle)
        vDelta_2 = (self.rise * (self.nT - 1)) - (webHeight + self.tf)
        flange_y = (self.w - self.tw) / 2
        front = -self.tT - self.tf
        outer = -self.tO - self.tw - flange_y

        coords = []
        if self.tp > 0:
            # Upper-Outer flange:
            coords.append(Vector([front, outer, -self.rise]))
            coords.append(Vector([-self.tT, outer, -self.rise]))
            coords.append(Vector([-self.tT, outer, 0]))
            coords.append(Vector([(self.run * (self.nT - 1)) - self.tT, outer,
                                  self.rise * (self.nT - 1)]))
            coords.append(Vector([self.run * self.nT, outer,
                                  self.rise * (self.nT - 1)]))
            coords.append(Vector([self.run * self.nT, outer,
                                  (self.rise * (self.nT - 1)) + self.tf]))
            coords.append(Vector([(self.run * (self.nT - 1)) - self.tT, outer,
                                  (self.rise * (self.nT - 1)) + self.tf]))
            coords.append(Vector([front, outer, self.tf - vDelta_1]))
            # Lower-Outer flange:
            coords.append(coords[0] + Vector([self.tf + webHeight, 0, 0]))
            coords.append(coords[1] + Vector([self.tf + webHeight, 0, 0]))
            coords.append(intersect_line_line(coords[9],
                                              coords[9] - Vector([0, 0, 1]),
                                              Vector([self.run, 0, -self.hT - self.tf]),
                                              Vector([self.run * 2, 0, self.rise - self.hT - self.tf]))[0])
            coords.append(Vector([(self.run * self.nT) - ((webHeight - self.hT) / math.tan(self.G.angle)),
                                  outer, vDelta_2]))
            coords.append(coords[4] - Vector([0, 0, self.tf + webHeight]))
            coords.append(coords[5] - Vector([0, 0, self.tf + webHeight]))
            coords.append(coords[11] + Vector([0, 0, self.tf]))
            coords.append(intersect_line_line(coords[8],
                                              coords[8] - Vector([0, 0, 1]),
                                              Vector([self.run, 0, -self.hT]),
                                              Vector([self.run * 2, 0, self.rise - self.hT]))[0])
            # Outer web:
            coords.append(coords[1] + Vector([0, flange_y, 0]))
            coords.append(coords[8] + Vector([0, flange_y, 0]))
            coords.append(coords[15] + Vector([0, flange_y, 0]))
            coords.append(coords[14] + Vector([0, flange_y, 0]))
            coords.append(coords[13] + Vector([0, flange_y, 0]))
            coords.append(coords[4] + Vector([0, flange_y, 0]))
            coords.append(coords[3] + Vector([0, flange_y, 0]))
            coords.append(coords[2] + Vector([0, flange_y, 0]))
            # Outer corner nodes:
            for i in [0, 7, 6, 5, 12, 11, 10, 9]:
                coords.append(coords[i] + Vector([0, flange_y + self.tw, 0]))

            self.G.Make_mesh(coords, self.faces4c, 'stringer')

            for i in range(16):
                coords[i] += Vector([0, -outer * 2, 0])
            for i in range(8):
                coords[i + 16] += Vector([0, (-outer - flange_y) * 2, 0])
            for i in coords:
                i += Vector([0, (self.tO * 2) + self.wT, 0])

            self.G.Make_mesh(coords, self.faces4c, 'stringer')

        return {'FINISHED'}

# from tread.py


class Treads:
    def __init__(self,G,typ,typ_t,run,w,h,d,r,toe,o,n,tk,sec,sp,sn,deg=4):
        self.G = G #General
        self.typ = typ #Stair type
        self.typ_t = typ_t #Tread type
        self.run = run #Stair run.  Degrees if self.typ == "id4"
        self.w=w #tread width.  Is outer radius if self.typ == "id4"
        self.h=h #tread height
        self.d=d #tread run.  Ignore for now if self.typ == "id4"
        self.r=r #tread rise
        self.t=toe #tread nosing
        self.o=o #tread side overhang.  Is inner radius if self.typ == "id4"
        self.n=n #number of treads
        self.tk=tk #thickness of tread metal
        self.sec=sec #metal sections for tread
        if sec != 1 and typ_t not in ["tId4", "tId5"]:
            self.sp=((d+toe)*(sp/100))/(sec-1) #spacing between sections (% of depth)
        elif typ_t in ["tId4", "tId5"]:
            self.sp=sp/100 #keep % value
        else:
            self.sp=0
        self.sn=sn #number of cross sections
        self.deg = deg #number of section per "slice".  Only applys if self.typ == "id4"
        self.tId2_faces = [[0,1,2,3],[0,3,4,5],[4,5,6,7],[6,7,8,9],[8,9,10,11],
                           [12,13,14,15],[12,15,16,17],[16,17,18,19],
                           [18,19,20,21],[20,21,22,23],[0,1,13,12],[1,2,14,13],
                           [2,3,15,14],[3,4,16,15],[4,7,19,16],[7,8,20,19],
                           [8,11,23,20],[11,10,22,23],[10,9,21,22],[9,6,18,21],
                           [6,5,17,18],[5,0,12,17]]
        self.out_faces = [[0,2,3,1],[0,2,10,8],[9,11,3,1],[9,11,10,8],
                          [2,6,7,3],[2,6,14,10],[11,15,7,3],[11,15,14,10],
                          [0,4,5,1],[0,4,12,8],[9,13,5,1],[9,13,12,8],
                          [4,6,7,5],[4,6,14,12],[13,15,14,12],[13,15,7,5]]
        self.Create()

    def Create(self):
        # Setup the coordinates:
        coords = []
        coords2 = []
        coords3 = []
        cross = 0
        cW = 0
        depth = 0
        offset = 0
        height = 0
        if self.typ in ["id1", "id2", "id3"]:
            if self.typ_t == "tId1":
                coords.append(Vector([-self.t,-self.o,0]))
                coords.append(Vector([self.d,-self.o,0]))
                coords.append(Vector([-self.t,self.w + self.o,0]))
                coords.append(Vector([self.d,self.w + self.o,0]))
                for i in range(4):
                    coords.append(coords[i]+Vector([0,0,-self.h]))

            elif self.typ_t == "tId2":
                depth = (self.d + self.t - (self.sec - 1) * self.sp) / self.sec
                inset = depth / 4
                tDepth = depth - self.t
                coords.append(Vector([-self.t, -self.o, -self.h]))                          #0
                coords.append(Vector([inset - self.t, -self.o, -self.h]))           #1
                coords.append(Vector([inset - self.t, -self.o, -self.h + self.tk])) #2
                coords.append(Vector([self.tk - self.t, -self.o, -self.h + self.tk]))       #3
                coords.append(Vector([self.tk - self.t, -self.o, -self.tk]))                #4
                coords.append(Vector([-self.t, -self.o, 0]))                                #5
                coords.append(Vector([tDepth, -self.o, 0]))                                 #6
                coords.append(Vector([tDepth - self.tk, -self.o, -self.tk]))                #7
                coords.append(Vector([tDepth - self.tk, -self.o, self.tk - self.h]))        #8
                coords.append(Vector([tDepth, -self.o, -self.h]))                           #9
                coords.append(Vector([tDepth - inset, -self.o, -self.h]))           #10
                coords.append(Vector([tDepth - inset, -self.o, -self.h + self.tk])) #11
                for i in range(12):
                    coords.append(coords[i] + Vector([0, self.w + (2 * self.o), 0]))

            elif self.typ_t in ["tId3", "tId4", "tId5"]:
                # Frame:
                coords.append(Vector([-self.t,-self.o,-self.h]))
                coords.append(Vector([self.d,-self.o,-self.h]))
                coords.append(Vector([-self.t,-self.o,0]))
                coords.append(Vector([self.d,-self.o,0]))
                for i in range(4):
                    if (i % 2) == 0:
                        coords.append(coords[i] + Vector([self.tk,self.tk,0]))
                    else:
                        coords.append(coords[i] + Vector([-self.tk,self.tk,0]))
                for i in range(4):
                    coords.append(coords[i] + Vector([0,self.w + self.o,0]))
                for i in range(4):
                    coords.append(coords[i + 4] + Vector([0,self.w + self.o - (2 * self.tk),0]))

                # Tread sections:
                if self.typ_t == "tId3":
                    offset = (self.tk * math.sqrt(2)) / 2
                    topset = self.h - offset
                    self.sp = ((self.d + self.t - (2 * self.tk)) - (offset * (self.sec) + topset)) / (self.sec + 1)
                    baseX = -self.t + self.sp + self.tk
                    coords2.append(Vector([baseX, self.tk - self.o, offset - self.h]))
                    coords2.append(Vector([baseX + offset, self.tk - self.o, -self.h]))
                    for i in range(2):
                        coords2.append(coords2[i] + Vector([topset, 0, topset]))
                    for i in range(4):
                        coords2.append(coords2[i] + Vector([0, (self.w + self.o) - (2 * self.tk), 0]))
                elif self.typ_t in ["tId4", "tId5"]:
                    offset = ((self.run + self.t) * self.sp) / (self.sec + 1)
                    topset = (((self.run + self.t) * (1 - self.sp)) - (2 * self.tk)) / self.sec
                    baseX = -self.t + self.tk + offset
                    baseY = self.w + self.o - 2 * self.tk
                    coords2.append(Vector([baseX, -self.o + self.tk, -self.h / 2]))
                    coords2.append(Vector([baseX + topset, -self.o + self.tk, -self.h / 2]))
                    coords2.append(Vector([baseX, -self.o + self.tk, 0]))
                    coords2.append(Vector([baseX + topset, -self.o + self.tk, 0]))
                    for i in range(4):
                        coords2.append(coords2[i] + Vector([0, baseY, 0]))

                # Tread cross-sections:
                if self.typ_t in ["tId3", "tId4"]:
                    cW = self.tk
                    cross = (self.w + (2 * self.o) - (self.sn + 2) * self.tk) / (self.sn + 1)
                else: # tId5
                    spacing = self.sp ** (1 / 4)
                    cross = ((2*self.o + self.w) * spacing) / (self.sn + 1)
                    cW = (-2*self.tk + (2*self.o + self.w) * (1 - spacing)) / self.sn
                    self.sp = topset
                    height = -self.h / 2
                baseY = -self.o + self.tk + cross
                coords3.append(Vector([-self.t + self.tk, baseY, -self.h]))
                coords3.append(Vector([self.d - self.tk, baseY, -self.h]))
                coords3.append(Vector([-self.t + self.tk, baseY, height]))
                coords3.append(Vector([self.d - self.tk, baseY, height]))
                for i in range(4):
                    coords3.append(coords3[i] + Vector([0, cW, 0]))

            # Make the treads:
            for i in range(self.n):
                if self.typ_t == "tId1":
                    self.G.Make_mesh(coords,self.G.faces,'treads')
                elif self.typ_t == "tId2":
                    temp = []
                    for j in coords:
                        temp.append(copy(j))
                    for j in range(self.sec):
                        self.G.Make_mesh(temp, self.tId2_faces, 'treads')
                        for k in temp:
                            k += Vector([depth + self.sp, 0, 0])
                elif self.typ_t in ["tId3", "tId4", "tId5"]:
                    self.G.Make_mesh(coords,self.out_faces,'treads')
                    temp = []
                    for j in coords2:
                        temp.append(copy(j))
                    for j in range(self.sec):
                        self.G.Make_mesh(temp,self.G.faces,'bars')
                        for k in temp:
                            k += Vector([offset + self.sp, 0, 0])
                    for j in coords2:
                        j += Vector([self.d, 0, self.r])
                    temp = []
                    for j in coords3:
                        temp.append(copy(j))
                    for j in range(self.sn):
                        self.G.Make_mesh(temp,self.G.faces,'crosses')
                        for k in temp:
                            k += Vector([0, cW + cross, 0])
                    for j in coords3:
                        j += Vector([self.d, 0, self.r])
                for j in coords:
                    j += Vector([self.d,0,self.r])
        # Circular staircase:
        elif self.typ in ["id4"]:
            start = [Vector([0, -self.o, 0]), Vector([0, -self.o, -self.h]),
                     Vector([0, -self.w, 0]), Vector([0, -self.w, -self.h])]
            self.d = math.radians(self.run) / self.n
            for i in range(self.n):
                coords = []
                # Base faces.  Should be able to append more sections:
                tId4_faces = [[0, 1, 3, 2]]
                t_inner = Matrix.Rotation((-self.t / self.o) + (self.d * i), 3, 'Z')
                coords.append((t_inner * start[0]) + Vector([0, 0, self.r * i]))
                coords.append((t_inner * start[1]) + Vector([0, 0, self.r * i]))
                t_outer = Matrix.Rotation((-self.t / self.w) + (self.d * i), 3, 'Z')
                coords.append((t_outer * start[2]) + Vector([0, 0, self.r * i]))
                coords.append((t_outer * start[3]) + Vector([0, 0, self.r * i]))
                k = 0
                for j in range(self.deg + 1):
                    k = (j * 4) + 4
                    tId4_faces.append([k, k - 4, k - 3, k + 1])
                    tId4_faces.append([k - 2, k - 1, k + 3, k + 2])
                    tId4_faces.append([k + 1, k - 3, k - 1, k + 3])
                    tId4_faces.append([k, k - 4, k - 2, k + 2])
                    rot = Matrix.Rotation(((self.d * j) / self.deg) + (self.d * i), 3, 'Z')
                    for v in start:
                        coords.append((rot * v) + Vector([0, 0, self.r * i]))
                tId4_faces.append([k, k + 1, k + 3, k + 2])
                self.G.Make_mesh(coords, tId4_faces, 'treads')
        return


class stairs(bpy.types.Operator):
    """Add stair objects"""
    bl_idname = "mesh.stairs"
    bl_label = "Add Stairs"
    bl_options = {'REGISTER', 'UNDO'}
    bl_description = "Add stairs"

    # Stair types for enum:
    id1 = ("id1", "Freestanding", "Generate a freestanding staircase.")
    id2 = ("id2", "Housed-Open", "Generate a housed-open staircase.")
    id3 = ("id3", "Box", "Generate a box staircase.")
    id4 = ("id4", "Circular", "Generate a circular or spiral staircase.")

    # Tread types for enum:
    tId1 = ("tId1", "Classic", "Generate wooden style treads")
    tId2 = ("tId2", "Basic Steel", "Generate common steel style treads")
    tId3 = ("tId3", "Bar 1", "Generate bar/slat steel treads")
    tId4 = ("tId4", "Bar 2", "Generate bar-grating steel treads")
    tId5 = ("tId5", "Bar 3", "Generate bar-support steel treads")

    # Stringer types for enum:
    sId1 = ("sId1", "Classic", "Generate a classic style stringer")
    sId2 = ("sId2", "I-Beam", "Generate a steel I-beam stringer")
    sId3 = ("sId3", "C-Beam", "Generate a C-channel style stringer")

    typ = EnumProperty(name = "Type",
                       description = "Type of staircase to generate",
                       items = [id1, id2, id3, id4])

    rise = FloatProperty(name = "Rise",
                         description = "Single tread rise",
                         min = 0.0, max = 1024.0,
                         default = 0.20)
    run = FloatProperty(name = "Run",
                        description = "Single tread run",
                        min = 0.0, max = 1024.0,
                        default = 0.30)

    #for circular
    rad1 = FloatProperty(name = "Inner Radius",
                         description = "Inner radius for circular staircase",
                         min = 0.0, max = 1024.0,
                         soft_max = 10.0,
                         default = 0.25)
    rad2 = FloatProperty(name = "Outer Radius",
                         description = "Outer radius for circular staircase",
                         min = 0.0, max = 1024.0,
                         soft_min = 0.015625, soft_max = 32.0,
                         default = 1.0)
    deg = FloatProperty(name = "Degrees",
                        description = "Number of degrees the stairway rotates",
                        min = 0.0, max = 92160.0, step = 5.0,
                        default = 450.0)
    center = BoolProperty(name = "Center Pillar",
                          description = "Generate a central pillar",
                          default = False)

    #for treads
    make_treads = BoolProperty(name = "Make Treads",
                              description = "Enable tread generation",
                              default = True)
    tread_w = FloatProperty(name = "Tread Width",
                            description = "Width of each generated tread",
                            min = 0.0001, max = 1024.0,
                            default = 1.2)
    tread_h = FloatProperty(name = "Tread Height",
                            description = "Height of each generated tread",
                            min = 0.0001, max = 1024.0,
                            default = 0.04)
    tread_t = FloatProperty(name = "Tread Toe",
                            description = "Toe (aka \"nosing\") of each generated tread",
                            min = 0.0, max = 10.0,
                            default = 0.03)
    tread_o = FloatProperty(name = "Tread Overhang",
                            description = "How much tread \"overhangs\" the sides",
                            min = 0.0, max = 1024.0,
                            default = 0.025)
    tread_n = IntProperty(name = "Number of Treads",
                          description = "How many treads to generate",
                          min = 1, max = 1024,
                          default = 10)
    typ_t = EnumProperty(name = "Tread Type",
                         description = "Type/style of treads to generate",
                         items = [tId1, tId2, tId3, tId4, tId5])
    tread_tk = FloatProperty(name = "Thickness",
                             description = "Thickness of the treads",
                             min = 0.0001, max = 10.0,
                             default = 0.02)
    tread_sec = IntProperty(name = "Sections",
                            description = "Number of sections to use for tread",
                            min = 1, max = 1024,
                            default = 5)
    tread_sp = IntProperty(name = "Spacing",
                           description = "Total spacing between tread sections as a percentage of total tread width",
                           min = 0, max = 80,
                           default = 5)
    tread_sn = IntProperty(name = "Crosses",
                           description = "Number of cross section supports",
                           min = 2, max = 1024,
                           default = 4)
    #special circular tread properties:
    tread_slc = IntProperty(name = "Slices",
                            description = "Number of slices each tread is composed of",
                            min = 1, max = 1024,
                            soft_max = 16,
                            default = 4)

    #for posts
    make_posts = BoolProperty(name = "Make Posts",
                              description = "Enable post generation",
                              default = True)
    post_d = FloatProperty(name = "Post Depth",
                           description = "Depth of generated posts",
                           min = 0.0001, max = 10.0,
                           default = 0.04)
    post_w = FloatProperty(name = "Post Width",
                           description = "Width of generated posts",
                           min = 0.0001, max = 10.0,
                           default = 0.04)
    post_n = IntProperty(name = "Number of Posts",
                         description = "Number of posts to generated",
                         min = 1, max = 1024,
                         default = 5)

    #for railings
    make_railings = BoolProperty(name = "Make Railings",
                                 description = "Generate railings",
                                 default = True)
    rail_w = FloatProperty(name = "Railings Width",
                           description = "Width of railings to generate",
                           min = 0.0001, max = 10.0,
                           default = 0.12)
    rail_t = FloatProperty(name = "Railings Thickness",
                           description = "Thickness of railings to generate",
                           min = 0.0001, max = 10.0,
                           default = 0.03)
    rail_h = FloatProperty(name = "Railings Height",
                           description = "Height of railings to generate",
                           min = 0.0001, max = 10.0,
                           default = 0.90)

    #for retainers
    make_retainers = BoolProperty(name = "Make Retainers",
                                  description = "Generate retainers",
                                  default = True)
    ret_w = FloatProperty(name = "Retainer Width",
                          description = "Width of generated retainers",
                          min = 0.0001, max = 10.0,
                          default = 0.01)
    ret_h = FloatProperty(name = "Retainer Height",
                          description = "Height of generated retainers",
                          min = 0.0001, max = 10.0,
                          default = 0.01)
    ret_n = IntProperty(name = "Number of Retainers",
                        description = "Number of retainers to generated",
                        min = 1, max = 1024,
                        default = 3)

    #for stringer
    make_stringer = BoolProperty(name = "Make Stringer",
                                 description = "Generate stair stringer",
                                 default = True)
    typ_s = EnumProperty(name = "Stringer Type",
                         description = "Type/style of stringer to generate",
                         items = [sId1, sId2, sId3])
    string_n = IntProperty(name = "Number of Stringers",
                           description = "Number of stringers to generate",
                           min = 1, max = 10,
                           default = 1)
    string_dis = BoolProperty(name = "Distributed",
                              description = "Use distributed stringers",
                              default = False)
    string_w = FloatProperty(name = "Stringer width",
                             description = "Width of stringer as a percentage of tread width",
                             min = 0.0001, max = 100.0,
                             default = 15.0)
    string_h = FloatProperty(name = "Stringer Height",
                             description = "Height of the stringer",
                             min = 0.0001, max = 100.0,
                             default = 0.3)
    string_tw = FloatProperty(name = "Web Thickness",
                              description = "Thickness of the beam's web as a percentage of width",
                              min = 0.0001, max = 100.0,
                              default = 25.0)
    string_tf = FloatProperty(name = "Flange Thickness",
                              description = "Thickness of the flange",
                              min = 0.0001, max = 100.0,
                              default = 0.05)
    string_tp = FloatProperty(name = "Flange Taper",
                              description = "Flange thickness taper as a percentage",
                              min = 0.0, max = 100.0,
                              default = 0.0)
    string_g = BoolProperty(name = "Floating",
                            description = "Cut bottom of strigner to be a \"floating\" section",
                            default = False)

    use_original = BoolProperty(name = "Use legacy method",
                                description = "Use the Blender 2.49 legacy method for stair generation",
                                default = True)
    rEnable = BoolProperty(name = "Right Details",
                           description = "Generate right side details (posts/rails/retainers)",
                           default = True)
    lEnable = BoolProperty(name = "Left Details",
                           description = "Generate left side details (posts/rails/retainers)",
                           default = True)

    # Draw the GUI:
    def draw(self, context):
        layout = self.layout
        box = layout.box()
        box.prop(self, 'typ')
        box = layout.box()
        box.prop(self, 'rise')
        if self.typ != "id4":
            box.prop(self, 'run')
        else:
            box.prop(self, 'deg')
            box.prop(self, 'rad1')
            box.prop(self, 'rad2')
            box.prop(self, 'center')
        if self.typ == "id1":
            box.prop(self, 'use_original')
            if not self.use_original:
                box.prop(self, 'rEnable')
                box.prop(self, 'lEnable')
        else:
            self.use_original = False
            box.prop(self, 'rEnable')
            box.prop(self, 'lEnable')

        # Treads
        box = layout.box()
        box.prop(self, 'make_treads')
        if self.make_treads:
            if not self.use_original and self.typ != "id4":
                box.prop(self, 'typ_t')
            else:
                self.typ_t = "tId1"
            if self.typ != "id4":
                box.prop(self, 'tread_w')
            box.prop(self, 'tread_h')
            box.prop(self, 'tread_t')
            if self.typ not in ["id2", "id4"]:
                box.prop(self, 'tread_o')
            else:
                self.tread_o = 0.0
            box.prop(self, 'tread_n')
            if self.typ_t != "tId1":
                box.prop(self, 'tread_tk')
                box.prop(self, 'tread_sec')
                if self.tread_sec > 1 and self.typ_t not in ["tId3", "tId4"]:
                    box.prop(self, 'tread_sp')
                if self.typ_t in ["tId3", "tId4", "tId5"]:
                    box.prop(self, 'tread_sn')
            elif self.typ == "id4":
                box.prop(self, "tread_slc")

        # Posts
        box = layout.box()
        box.prop(self, 'make_posts')
        if self.make_posts:
            box.prop(self, 'post_d')
            box.prop(self, 'post_w')
            box.prop(self, 'post_n')

        # Railings
        box = layout.box()
        box.prop(self, 'make_railings')
        if self.make_railings:
            box.prop(self, 'rail_w')
            box.prop(self, 'rail_t')
            box.prop(self, 'rail_h')

        # Retainers
        box = layout.box()
        box.prop(self, 'make_retainers')
        if self.make_retainers:
            box.prop(self, 'ret_w')
            box.prop(self, 'ret_h')
            box.prop(self, 'ret_n')

        # Stringers
        box = layout.box()
        if self.typ != "id2":
            box.prop(self, 'make_stringer')
        else:
            self.make_stringer = True
        if self.make_stringer:
            if not self.use_original:
                box.prop(self, 'typ_s')
            else:
                self.typ_s = "sId1"
            box.prop(self, 'string_w')
            if self.typ == "id1":
                if self.typ_s == "sId1" and not self.use_original:
                    box.prop(self, 'string_n')
                    box.prop(self, 'string_dis')
                elif self.typ_s in ["sId2", "sId3"]:
                    box.prop(self, 'string_n')
                    box.prop(self, 'string_dis')
                    box.prop(self, 'string_h')
                    box.prop(self, 'string_tw')
                    box.prop(self, 'string_tf')
                    box.prop(self, 'string_tp')
                    box.prop(self, 'string_g')
            elif self.typ == "id2":
                if self.typ_s in ["sId2", "sId3"]:
                    box.prop(self, 'string_tw')
                    box.prop(self, 'string_tf')

        # Tread support:
##        if self.make_stringer and typ_s in ["sId2", "sId3"]:

    def execute(self, context):
        global G
        global typ
        global typ_s
        global typ_t
        global rise
        global run
        typ = self.typ
        typ_s = self.typ_s
        typ_t = self.typ_t
        rise = self.rise
        run = self.run
        G=General(rise,run,self.tread_n)
        if self.make_treads:
            if typ != "id4":
                Treads(G,
                       typ,
                       typ_t,
                       run,
                       self.tread_w,
                       self.tread_h,
                       self.run,
                       self.rise,
                       self.tread_t,
                       self.tread_o,
                       self.tread_n,
                       self.tread_tk,
                       self.tread_sec,
                       self.tread_sp,
                       self.tread_sn)
            else:
                Treads(G,
                       typ,
                       typ_t,
                       self.deg,
                       self.rad2,
                       self.tread_h,
                       self.run,
                       self.rise,
                       self.tread_t,
                       self.rad1,
                       self.tread_n,
                       self.tread_tk,
                       self.tread_sec,
                       self.tread_sp,
                       self.tread_sn,
                       self.tread_slc)
        if self.make_posts and (self.rEnable or self.lEnable):
            Posts(G,
                  rise,
                  run,
                  self.post_d,
                  self.post_w,
                  self.tread_w,
                  self.post_n,
                  self.rail_h,
                  self.rail_t,
                  self.rEnable,
                  self.lEnable)
        if self.make_railings and (self.rEnable or self.lEnable):
            Rails(G,
                  self.rail_w,
                  self.rail_t,
                  self.rail_h,
                  self.tread_t,
                  self.post_w,
                  self.post_d,
                  self.tread_w,
                  self.rEnable,
                  self.lEnable)
        if self.make_retainers and (self.rEnable or self.lEnable):
            Retainers(G,
                      self.ret_w,
                      self.ret_h,
                      self.post_w,
                      self.tread_w,
                      self.rail_h,
                      self.ret_n,
                      self.rEnable,
                      self.lEnable)
        if self.make_stringer:
            if typ == "id1" and self.use_original:
                Stringer(G,
                         typ,
                         typ_s,
                         rise,
                         run,
                         self.string_w,
                         self.string_h,
                         self.tread_n,
                         self.tread_h,
                         self.tread_w,
                         self.tread_t,
                         self.tread_o,
                         self.string_tw,
                         self.string_tf,
                         self.string_tp,
                         not self.string_g)
            elif typ == "id3":
                Stringer(G,
                         typ,
                         typ_s,
                         rise,
                         run,
                         100,
                         self.string_h,
                         self.tread_n,
                         self.tread_h,
                         self.tread_w,
                         self.tread_t,
                         self.tread_o,
                         self.string_tw,
                         self.string_tf,
                         self.string_tp,
                         not self.string_g,
                         1, False, False)
            elif typ == "id4":
                Stringer(G,
                         typ,
                         typ_s,
                         rise,
                         self.deg,
                         self.string_w,
                         self.string_h,
                         self.tread_n,
                         self.tread_h,
                         self.rad2 - self.rad1,
                         self.tread_t,
                         self.rad1,
                         self.string_tw,
                         self.string_tf,
                         self.string_tp,
                         not self.string_g,
                         self.string_n,
                         self.string_dis,
                         self.use_original,
                         self.tread_slc)
            else:
                Stringer(G,
                         typ,
                         typ_s,
                         rise,
                         run,
                         self.string_w,
                         self.string_h,
                         self.tread_n,
                         self.tread_h,
                         self.tread_w,
                         self.tread_t,
                         self.tread_o,
                         self.string_tw,
                         self.string_tf,
                         self.string_tp,
                         not self.string_g,
                         self.string_n,
                         self.string_dis,
                         self.use_original)
        return {'FINISHED'}
