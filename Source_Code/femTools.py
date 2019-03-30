# ***************************************************************************
# *                                                                         *
# *   Copyright (c) 2019 - Harry van Langen <hvlanalysis@icloud.com>        *
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU Lesser General Public License (LGPL)    *
# *   as published by the Free Software Foundation; either version 2 of     *
# *   the License, or (at your option) any later version.                   *
# *   for detail see the LICENCE text file.                                 *
# *                                                                         *
# *   This program is distributed in the hope that it will be useful,       *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU Library General Public License for more details.                  *
# *                                                                         *
# *   You should have received a copy of the GNU Library General Public     *
# *   License along with this program; if not, write to the Free Software   *
# *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
# *   USA                                                                   *
# *                                                                         *
# ***************************************************************************

import FreeCAD as App
import FreeCADGui as Gui
import numpy as np
import ObjectsFem
import FemGui
from femmesh import meshtools as mt
from feminout import importToolsFem as itf
import femsolver.calculix.writer as ccxw
import ObjectsFem as objectsfem
from femobjects import _FemMaterial
import femtools.ccxtools as tools
import Part as Part
import sys
from femtools.femutils import get_several_member as gsm
from femsolver.writerbase import FemInputWriter as iw
import DraftVecUtils as DVU
import scipy.linalg
import math
import Fem
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import os

np.set_printoptions(precision=5, linewidth=300)

def setUpAnalysis():

    doc = App.ActiveDocument

    mesh = doc.getObject("FEMMeshGmsh").FemMesh
    if mesh == None:
        print("No Gmsh object. Please create one first")
        raise SystemExit ()

    analysis = doc.getObject("Analysis")
    if analysis == None:
        print("No Analysis object. Please create one first")
        raise SystemExit ()

    # purge result objects
    for obj in App.ActiveDocument.Objects:
        name = obj.Name[:11]
        if name in ['MechanicalR', 'Result_Mesh']:
            doc.removeObject(obj.Name)

    doc.recompute()

    return doc, mesh, analysis

def setUpInput(doc, mesh, analysis):

    analysis = doc.getObject("Analysis")

    # create connextivity array elnodes for mapping local node number -> global node number
    elnodes = np.array([mesh.getElementNodes(el) for el in mesh.Volumes]) # elnodes[elementIndex] = [node1,...,Node10]
    elo = dict(zip(mesh.Volumes, range(len(mesh.Volumes)))) # elo : {elementNumber : elementIndex}

    # create nodal coordinate array nocoord for node number -> (x,y,z)
    nocoord=np.asarray(mesh.Nodes.values()) # nocoord[nodeIndex] = [x-coord, y-coord, z-coord]

    # create element material array: materialbyElement maps element number -> E, nu
    materials_lin = gsm(analysis, 'Fem::Material')
    fiwc = iw(
        analysis, doc.CalculiXccxTools, doc.FEMMeshGmsh, materials_lin,
        None, None, None, None, None, None, None, None,
        None, None, None, None, None, None, None, None, None
    )
    fiwc.get_material_elements()
    materialbyElement = []
    po_keys=[] # parentobject[el] = parent material object for element el
    po_values=[]
    counter=0
    for object in fiwc.material_objects:
        E = float(App.Units.Quantity(object['Object'].Material['YoungsModulus']).getValueAs('MPa'))
        Nu = float(object['Object'].Material['PoissonRatio'])
        for el in object['FEMElements']:
            po_keys.append(el)
            po_values.append(object['Object'].Name)
            counter += 1
            materialbyElement.append([el, E, Nu]) # materialbyElement[elementIndex] = [elementNumber, E, Nu]
    parentobject=dict(zip(po_keys,po_values))

    # determine elements connected to a node using FC API
    fet=mt.get_femelement_table(mesh) # fet is dictionary: { elementid : [ nodeid, nodeid, ... , nodeid ] }
    net=mt.get_femnodes_ele_table(mesh.Nodes, fet) # net is dictionary: {nodeID : [[eleID, NodePosition], [], ...], nodeID : [[], [], ...], ...}

    # set up interface element connectivity
    nodecount=len(nocoord)
    twins={} # twins : {oldnode : [newnode, face number]}
    interface_elements=[]
    for obj in App.ActiveDocument.Objects:
        if obj.Name == "BooleanFragments":
            num_BF_els=len(App.ActiveDocument.BooleanFragments.Objects)
            num_Mat_obs=len(fiwc.material_objects)
            if num_BF_els != num_Mat_obs:
                print("Each BooleanFragment element needs its own material object")
                raise SystemExit()
            shapecontacts=[]
            tk=[]
            tv=[]
            contactnum=0 # number of contact faces between BooleanFragment elements
            for index, obj1 in enumerate(obj.Objects):
                for obj1face in obj1.Shape.Faces:
                    el1faces=np.asarray(mesh.getFacesByFace(obj1face))
                    for obj2 in obj.Objects[index+1:]:
                        for obj2face in obj2.Shape.Faces:
                            el2faces = np.asarray(mesh.getFacesByFace(obj2face))
                            contact=np.intersect1d(el1faces,el2faces) # all volume element contact faces between obj1 and obj2
                            if contact != []:
                                contactnum+=1
                                npflist=[]
                                shapecontacts.append(contact) # all volume element contact faces
                                for contactface in contact:
                                    npflist.append(mesh.getElementNodes(contactface)) # all nodes at the contact between obj1 and obj2
                                # flatten npflist and remove duplicate nodes
                                cn = list(set([node for npf in npflist for node in
                                             npf]))
                                cn.sort()
                                for node in cn:
                                    if node not in tk:
                                        nodecount+=1
                                        tk.append(int(node)) # add an existing contact node to interface dict
                                        tv.append([nodecount, contactnum]) # add a new contact node (nodecount) to interface dict and the contact (contactnum) it is a member of
                                        nocoord = np.append(nocoord, [nocoord[node-1]], axis=0) # add coordinates for new nodes

            twins = dict(zip(tk, tv)) # twins : {oldnode : [newnode, face number]}
            print("\nInterface twins: {}".format(twins))

            for facecontacts in shapecontacts:
                for face in facecontacts:
                    nodeset1=list(mesh.getElementNodes(face))
                    nodeset2=[]
                    for node in mesh.getElementNodes(face):
                        nodeset2.append(twins[node][0])
                    interface_elements.append(nodeset1+nodeset2) # interface_elements[index] = [oldnode1,..,oldnode6, newnode1,..,newnode6]

            interface_elements=np.asarray(interface_elements) # add interface elements for the face between obj1 and obj2

    print("number of interface elements: {}".format(len(interface_elements)))
    for ie in interface_elements:
        print(ie)

    # reconnect volume elements to new interface nodes
    if interface_elements != []:
        membership = [""] * contactnum
        shiftelements=[] # reconnected elements
        for node in twins:
            if membership[twins[node][1]-1] == "":
                membership[twins[node][1]-1] = parentobject[net[node][0][0]] # elements that will stay with old nodes are in membership material object
            for element in net[node]:
                elnum = element[0] # element number
                if parentobject[elnum]!=membership[twins[node][1]-1]: # reconnect element to new nodes
                    nonum = int(math.log(element[1],
                                         2))  # local node number from binary node number net[node][1]
                    shiftelements.append(elo[elnum])
                    elnodes[elo[elnum]][nonum] = twins[node][0] # connect element to new nodes

    noce=np.zeros((len(nocoord)), dtype=np.int16)
    for inno in range (len(nocoord)):
        i, j = np.where(elnodes == inno+1)
        noce[inno]=len(i)

    # create boundary condition array dispfaces
    dispfaces=[]
    for obj in App.ActiveDocument.Objects:
        if obj.isDerivedFrom('Fem::ConstraintDisplacement'):
            bcnodes= []
            bctype=[obj.xFree,obj.yFree,obj.zFree]
            bcvalue=[obj.xDisplacement,obj.yDisplacement,obj.zDisplacement]
            for part, boundaries in obj.References:
                for boundary in boundaries:
                    ref = part.Shape.getElement(boundary)
                    if type(ref) == Part.Vertex:
                        bc=mesh.getNodesByVertex(ref)
                        for bcn in bc: bcnodes.append(bcn)
                    elif type(ref) == Part.Edge:
                        bc=mesh.getNodesByEdge(ref)
                        for bcn in bc: bcnodes.append(bcn)
                    elif type(ref) == Part.Face:
                        bc=mesh.getNodesByFace(ref)
                        for bcn in bc:
                            if bcn not in twins:
                                bcnodes.append(bcn)
                    else:
                        print("No Boundaries Found")
                    edge=list(set(bc) & set(twins))
                    interior=list(set(bc) - set(twins))
                    #print("\nset(bc) {}, \nset(twins) {}, \nedge of boundary: {}, \ninterior of boundary: {}".format(bc, twins, edge, interior))
                    for node in edge:
                        for elem in net[node]:
                            elnum=elo[elem[0]]
                            if list(set(elnodes[elnum]) & set(interior)) != []:
                                if elnum in shiftelements:
                                    bcnodes.append(twins[node][0])
                                else:
                                    bcnodes.append(node)

            bcnodes = list(dict.fromkeys(bcnodes)) #remove duplicates in bcnodes
            if bcnodes != []: dispfaces.append([bcnodes,bctype,bcvalue])

    # create loaded faces and their global node numbers
    loadfaces_keys = []
    loadfaces_values = []
    for obj in App.ActiveDocument.Objects:
        if obj.isDerivedFrom('Fem::ConstraintPressure'):
            if obj.Reversed:
                sign=1
            else:
                sign=-1
            for part, boundaries in obj.References:
                for boundary in boundaries:
                    ref = part.Shape.getElement(boundary)
                    if type(ref)==Part.Face:
                        for faceID in mesh.getFacesByFace(ref):
                            loadfaces_keys.append(faceID)
                            loadfaces_values.append([list(mesh.getElementNodes(faceID)),sign*obj.Pressure])
                    else:
                        print("No Faces with Pressure Loads")
    loadfaces=dict(zip(loadfaces_keys,loadfaces_values))


    # re-order element nodes
    for el in elnodes:
        temp  = el[1]
        el[1] = el[2]
        el[2] = temp
        temp  = el[4]
        el[4] = el[6]
        el[6] = temp
        temp  = el[8]
        el[8] = el[9]
        el[9] = temp

    return elnodes, nocoord, dispfaces, loadfaces, materialbyElement, interface_elements, noce

# shape functions for a 4-node tetrahedron - only used for stress interpolation
def shape4tet(xi, et, ze, xl):
    shp = np.zeros((4), dtype=np.float64)

    # shape functions
    shp[0] = 1.0 - xi - et - ze
    shp[1] = xi
    shp[2] = et
    shp[3] = ze

    return shp

# shape functions and their derivatives for a 10-node tetrahedron
def shape10tet(xi, et, ze, xl):
    shp = np.zeros((10), dtype=np.float64)
    dshp = np.zeros((3, 10), dtype=np.float64)
    bmat = np.zeros((6, 30), dtype=np.float64)

    # shape functions - source: Calculix, G Dhondt
    a = 1.0 - xi - et - ze
    shp[0] = (2.0 * a - 1.0) * a
    shp[1] = xi * (2.0 * xi - 1.0)
    shp[2] = et * (2.0 * et - 1.0)
    shp[3] = ze * (2.0 * ze - 1.0)
    shp[4] = 4.0 * xi * a
    shp[5] = 4.0 * xi * et
    shp[6] = 4.0 * et * a
    shp[7] = 4.0 * ze * a
    shp[8] = 4.0 * xi * ze
    shp[9] = 4.0 * et * ze

    # local derivatives of the shape functions: xi-derivative - source: Calculix, G Dhondt
    dshp[0][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
    dshp[0][1] = 4.0 * xi - 1.0
    dshp[0][2] = 0.0
    dshp[0][3] = 0.0
    dshp[0][4] = 4.0 * (1.0 - 2.0 * xi - et - ze)
    dshp[0][5] = 4.0 * et
    dshp[0][6] = -4.0 * et
    dshp[0][7] = -4.0 * ze
    dshp[0][8] = 4.0 * ze
    dshp[0][9] = 0.0

    # local derivatives of the shape functions: eta-derivative - source: Calculix, G Dhondt
    dshp[1][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
    dshp[1][1] = 0.0
    dshp[1][2] = 4.0 * et - 1.0
    dshp[1][3] = 0.0
    dshp[1][4] = -4.0 * xi
    dshp[1][5] = 4.0 * xi
    dshp[1][6] = 4.0 * (1.0 - xi - 2.0 * et - ze)
    dshp[1][7] = -4.0 * ze
    dshp[1][8] = 0.0
    dshp[1][9] = 4.0 * ze

    # local derivatives of the shape functions: zeta-derivative - source: Calculix, G Dhondt
    dshp[2][0] = 1.0 - 4.0 * (1.0 - xi - et - ze)
    dshp[2][1] = 0.0
    dshp[2][2] = 0.0
    dshp[2][3] = 4.0 * ze - 1.0
    dshp[2][4] = -4.0 * xi
    dshp[2][5] = 0.0
    dshp[2][6] = -4.0 * et
    dshp[2][7] = 4.0 * (1.0 - xi - et - 2.0 * ze)
    dshp[2][8] = 4.0 * xi
    dshp[2][9] = 4.0 * et

    xs = np.dot(xl, dshp.T) # local derivative of the global coordinates

    xsj = np.linalg.det(xs) # Jacobian

    xsi = np.linalg.inv(xs) # global derivative of the local coordinates

    dshp = np.dot(xsi.T, dshp) # global derivatives of the shape functions

    # computation of the strain interpolation matrix bmat
    for i in range(10):
        i3=3*i
        d00 = dshp [0][i]
        d10 = dshp [1][i]
        d20 = dshp [2][i]
        bmat[0][i3] = d00
        bmat[1][i3+1] = d10
        bmat[2][i3+2] = d20
        bmat[3][i3] = d10
        bmat[3][i3+1] = d00
        bmat[4][i3] = d20
        bmat[4][i3+2] = d00
        bmat[5][i3+1] = d20
        bmat[5][i3+2] = d10

    return xsj, shp, dshp, bmat

# shape functions and their derivatives for a 6-node triangular interface element
def shape6tri(xi, et, xl):
    shp = np.zeros((6), dtype=np.float64)
    dshp = np.zeros((2, 6), dtype=np.float64)
    bmat = np.zeros((3, 36), dtype=np.float64)

    # shape functions
    shp[0] = (1.0-xi-et)*(1.0-2.0*xi-2.0*et)
    shp[1] = xi*(2.0 *xi-1.0)
    shp[2] = et*(2.0*et-1.0)
    shp[3] = 4.0*xi*(1.0-xi-et)
    shp[4] = 4.0*xi*et
    shp[5] = 4.0*et*(1-xi-et)

    # local derivatives of the shape functions: xi-derivative
    dshp[0][0] = -3.0+4.0*et+4.0*xi
    dshp[0][1] = -1.0+4.0*xi
    dshp[0][2] = 0.0
    dshp[0][3] = -4.0*(-1.0+et+2.0*xi)
    dshp[0][4] = 4.0*et
    dshp[0][5] = -4.0*et

    # local derivatives of the shape functions: eta-derivative
    dshp[1][0] = -3.0+4.0*et+4.0*xi
    dshp[1][1] = 0.0
    dshp[1][2] = -1.0+4.0*et
    dshp[1][3] =-4.0*xi
    dshp[1][4] = 4.0*xi
    dshp[1][5] =-4.0*(-1.0+2.0*et+xi)

    xs = np.dot(dshp, xl.T) # xs = [ [[dx/dxi],[dy/dxi],[dz/dxi]] , [[dx/det],[dy/det],[dz/det]] ]

    xp = np.cross(xs[0],xs[1]) # vector normal to surface

    xsj = np.linalg.norm(xp) # Jacobian

    xx = xs[0]/np.linalg.norm(xs[0]) # unit vector in xi direction
    xp /= xsj # unit vector normal to surface
    xt = np.cross(xp,xx) # unit vector tangential to surface and normal to xx

    # computation of the "strain" interpolation matrix bmat
    for i in range(6):
        ia=3*i
        ib=ia+18
        ni = shp [i]
        bmat[0][ia] = ni
        bmat[1][ia+1] = ni
        bmat[2][ia+2] = ni
        bmat[0][ib] = -ni
        bmat[1][ib+1] = -ni
        bmat[2][ib+2] = -ni

    return xsj, shp, bmat, xx, xt, xp

# linear-elastic material stiffness matrix
def hooke(element, materialbyElement):
    dmat = np.zeros((6, 6), dtype=np.float64)
    e = materialbyElement[element - 1][1] # Young's Modulus
    nu = materialbyElement[element - 1][2] # Poisson's Ratio
    dm = e * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu)
    od = nu / (1.0 - nu)
    sd = 0.5* (1.0 - 2.0 * nu) / (1.0 - nu)
    dmat[0][0] = dmat[1][1] = dmat[2][2] = 1.0
    dmat[3][3] = dmat[4][4] = dmat[5][5] = sd
    dmat[0][1] = dmat[0][2] = dmat[1][2] = od
    dmat[1][0] = dmat[2][0] = dmat[2][1] = od
    dmat *= dm
    return dmat

# Gaussian integration points and weights
def gaussPoints():
    # Gaussian integration points and weights for 10-noded tetrahedron
    gp10 = np.array([[0.138196601125011, 0.138196601125011, 0.138196601125011,
                     0.041666666666667],
                    [0.585410196624968, 0.138196601125011, 0.138196601125011,
                     0.041666666666667],
                    [0.138196601125011, 0.585410196624968, 0.138196601125011,
                     0.041666666666667],
                    [0.138196601125011, 0.138196601125011, 0.585410196624968,
                     0.041666666666667]])
    # Gaussian integration points and weights for 6-noded triangle
    gp6 = np.array([[0.445948490915965,0.445948490915965,
                      0.111690794839005],
                     [1.0-2.0*0.445948490915965,0.445948490915965,
                      0.111690794839005],
                     [0.445948490915965,1.0-2.0*0.445948490915965,
                      0.111690794839005],
                     [0.091576213509771,0.091576213509771,
                      0.054975871827661],
                     [1.0-2.0*0.091576213509771,0.091576213509771,
                      0.054975871827661],
                     [0.091576213509771,1.0-2.0*0.091576213509771,
                      0.054975871827661]])
    return gp10, gp6

# Nodal point locations
def nodalPoints():
    # Nodal point locations for a 10-noded tetrahedron
    np10 = np.array([[0.0, 0.0, 0.0],
                     [1.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0],
                     [0.0, 0.0, 1.0],
                     [0.5, 0.0, 0.0],
                     [0.5, 0.5, 0.0],
                     [0.0, 0.5, 0.0],
                     [0.0, 0.0, 0.5],
                     [0.5, 0.0, 0.5],
                     [0.0, 0.5, 0.5]])
    # Nodal point locations for 6-noded triangle + Newton Cotes integration weights
    np6 = np.array([[0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.5, 0.0, 0.166666666666666],
                    [0.5, 0.5, 0.166666666666666],
                    [0.0, 0.5, 0.166666666666666]])
    return np10, np6


# caculate the global stiffness matrix and load vector
def calcGSM(elnodes, nocoord, materialbyElement, loadfaces, interface_elements, grav, kn, ks):

    gp10, gp6 = gaussPoints()
    np10, np6 = nodalPoints()
    nn = len(nocoord[:, 0])
    # global stiffness matrix
    gsm = np.zeros((3*nn, 3*nn), dtype=np.float64)
    # global load vector
    glv = np.zeros((3*nn), dtype=np.float64)

    #   calculate element load vectors for pressure and add to global vector
    for face in loadfaces:
        pressure = loadfaces[face][1]
        xl = np.array([nocoord[nd-1] for nd in loadfaces[face][0]]).T

        # integrate element load vector
        for ip in gp6:
            xi = ip[0]
            et = ip[1]
            xsj, shp, bmat, xx, xt, xp = shape6tri(xi, et, xl)

            nl = 0
            for nd in loadfaces[face][0]:
                iglob = nd - 1
                iglob3 = 3 * iglob
                for k in range(3):
                    load=shp [nl] * pressure * xp [k] * abs(xsj) * ip[2]
                    glv [iglob3+k] += load
                nl+=1

    # for each volume element calculate the element stiffness matrix
    # and gravity load vector and add to global matrix and vector
    element = 0
    for nodes in elnodes:
        V=0.0
        element += 1
        esm = np.zeros((30, 30), dtype=np.float64)
        gamma = np.zeros((30), dtype=np.float64)
        # stress-strain matrix dmat
        dmat = hooke(element, materialbyElement)
        # set up nodal values for this element
        xl = np.array([nocoord[nd-1] for nd in nodes]).T
        # integrate element matrix
        for ip in gp10:
            xi = ip[0]
            et = ip[1]
            ze = ip[2]
            xsj, shp, dshp, bmat = shape10tet(xi, et, ze, xl)
            esm += np.dot(bmat.T, np.dot(dmat,bmat)) * ip[3] * abs(xsj)
            gamma [2::3] += grav * shp * ip[3] * abs(xsj)
            V+=xsj*ip[3] # Element volume - not used
        # add element matrix to global stiffness matrix and element gravity
        # load vector to global load vector
        for i in range(10):
            iglob = nodes[i]-1
            iglob3 = 3*iglob
            i3 = 3*i
            glv[iglob3+2] += gamma[i3+2]
            for j in range(10):
                jglob = nodes[j]-1
                jglob3 = 3*jglob
                j3 = j*3
                for k in range (3):
                    for l in range (3):
                        gsm[iglob3+k, jglob3+l] += esm[i3+k, j3+l]

    # interface stiffness value
    kmax=0.01*np.amax(np.diag(gsm))

    # For each interface element calculate the element matrix and
    # add to global stiffness matrix
    for nodes in interface_elements:
        xl = np.array([nocoord[nd-1] for nd in nodes[:6]]).T
        esm = np.zeros((36, 36), dtype=np.float64)
        dmatloc = np.diag([kn*kmax,ks*kmax,ks*kmax])
        # integrate element matrix (np6: Newton Cotes, gp6: Gauss)
        for ip in np6:
            xi = ip[0]
            et = ip[1]
            xsj, shp, bmat, xx, xt, xp = shape6tri(xi, et, xl)
            T=np.array([xp, xx, xt])
            dmatglob=np.dot(T.T,np.dot(dmatloc,T))
            esm += np.dot(bmat.T,np.dot(dmatglob,bmat)) * abs(xsj) * ip[2]
        # add Element matrix to global stiffness matrix
        for i in range(12):
            iglob = nodes[i]-1
            iglob3 = 3*iglob
            i3 = 3*i
            for j in range(12):
                jglob = nodes[j]-1
                jglob3 = 3*jglob
                j3 = j*3
                for k in range (3):
                    for l in range (3):
                        gsm[iglob3+k, jglob3+l] += esm[i3+k, j3+l]

    loadsumx=0.0
    loadsumy=0.0
    loadsumz=0.0
    for node in range(nn):
        dof = 3*node
        loadsumx+=glv[dof]
        loadsumy+=glv[dof+1]
        loadsumz+=glv[dof+2]
    print("\nsumFx {} sumFy {} sumFz {}".format(loadsumx,loadsumy, loadsumz))

    return gsm, glv, kmax

# calculate load-deflection curve
def calcDisp (elnodes, nocoord, dispfaces, materialbyElement, interface_elements, kmax, gsm, glv, nstep, iterat_max,
                                                     error_max, relax, scale_re,
                                                     scale_up, scale_dn,
                                                     sig_yield, shr_yield, kn,
                                                     ks, out_disp):
    import time

    ndof=len(glv)
    nelem = len(elnodes)
    ninter = len(interface_elements)


    if np.min(np.diag(gsm)) <= 0.0:
        print("non poitive definite matrix - check input")
        for i in range(ndof):
            if gsm[i, i] == 0.0: print(
                "DOF: {}; Coord: {} not attached".format(i, nocoord[i / 3]))
        raise SystemExit()

    # glv will be impacted by non-zero prescribed displacements, so make a copy
    # to preserve the external load vector for use in the iterative scheme
    qex=np.copy(glv)
    qnorm = np.linalg.norm(qex)
    if qnorm < 1.0: qnorm = 1.0

    # modify the global stiffness matrix and load vector for displacement BC
    gsm, glv, fixdof = bcGSM(gsm, glv, dispfaces)

    # Cholesky decomposition of the global stiffness matrix and elastic solution
    # TODO: Apply reverse Cuthill McKee and banded Cholesky to speed things up
    # TODO: Experiment with Intel Distributin for Python (Math Kernal Library) to optimize speed
    t0 = time.time()
    L=scipy.linalg.cho_factor(gsm,True,True,False)
    t1 = time.time()
    ue=scipy.linalg.cho_solve(L,glv,True,False)
    t2 = time.time()
    print("Cholesky Decomposition: {} s, Elastic Solution: {} s".format(t1-t0,t2-t1))

    # initiate analysis
    dl0=1.0/nstep
    dl=dl0
    du=dl*ue

    sig = np.array([np.zeros((24*nelem), dtype=np.float64)]) # stress in Tet10
    trac = np.array([np.zeros((18*ninter), dtype=np.float64)]) # contact stress in Tri6
    disp = np.array([np.zeros((ndof), dtype=np.float64)]) # displacement results
    lbd = np.zeros((1), dtype=np.float64) # load level

    step = -1
    cnt = True

    while (cnt == True):
        for istep in (range(nstep)):
            step += 1
            restart=0
            print("\nStep: {}".format(step))
            # Riks control vector
            a=du
            # lbd = load level
            lbd=np.append(lbd, lbd[step]+dl)
            # update stresses
            sig_update=update_stress(elnodes, nocoord, materialbyElement, sig_yield, sig[step], du)
            sig=np.append(sig, sig_update, axis=0)
            # update interface stresses
            trac_update = update_traction(interface_elements, nocoord, trac[step], du,
                                             kmax, kn, ks, shr_yield)
            trac=np.append(trac, trac_update, axis = 0)
            # calculate internal load vector
            qin=update_load(elnodes, interface_elements, nocoord, sig[step+1], trac[step+1])
            # calculate residual load vector
            r=fixdof * (lbd[step+1] * qex - qin)
            rnorm=np.linalg.norm(r)
            # out-of-balance error
            error=rnorm/qnorm
            iterat=0
            print("Iteration: {}, Error: {}".format(iterat, error))

            while error>error_max:
                iterat+=1
                # displacement corrrection
                due = scipy.linalg.cho_solve(L, relax*r, True, False)
                # Riks control correction to load level increment
                dl=-np.dot(a,due)/np.dot(a,ue)
                lbd[step+1]+=dl
                # Riks control correction to displacement increment
                du+=due+dl*ue
                # update stresses
                sig[step+1]=update_stress(elnodes, nocoord, materialbyElement, sig_yield, sig[step], du)
                # update interface stresses
                trac[step + 1] = update_traction(interface_elements, nocoord, trac[step], du, kmax, kn, ks, shr_yield)
                # calculate internal load vector
                qin=update_load(elnodes, interface_elements, nocoord, sig[step+1], trac[step+1])
                # calculate out of balance error
                r = fixdof * (lbd[step + 1] * qex - qin)
                rnorm = np.linalg.norm(r)
                error = rnorm / qnorm
                print("Iteration: {}, Error: {}".format(iterat,error))
                if iterat>iterat_max:
                    # scale down
                    if restart == 4: raise SystemExit()
                    restart+=1
                    if step>0:
                        dl = (lbd[step ] - lbd[step-1])/scale_re/restart
                        du = (disp[step ] - disp[step-1])/scale_re/restart
                    else:
                        # for first step only
                        dl=dl0/scale_re/restart
                        du = dl*ue/scale_re/restart
                    lbd[step+1] = lbd[step] + dl
                    sig[step + 1] = update_stress(elnodes, nocoord,
                                                  materialbyElement, sig_yield, sig[step], du)
                    trac[step + 1] = update_traction(interface_elements,
                                                     nocoord, trac[step], du,
                                                     kmax, kn, ks, shr_yield)
                    qin = update_load(elnodes, interface_elements, nocoord, sig[step + 1], trac[step+1])
                    r = fixdof * (lbd[step + 1] * qex - qin)
                    rnorm = np.linalg.norm(r)
                    error = rnorm / qnorm
                    iterat=0
            # update results at end of converged load step
            disp=np.append(disp,[disp[step]+du],axis=0)
            dl=lbd[step+1]-lbd[step]
            if iterat>10:
                # scale down
                dl/=scale_dn
                du/=scale_dn
            if iterat<5:
                #scale up
                dl*=scale_up
                du*=scale_up

        # maximum displacement increment for plotting load-displacement curve
        un=[]
        for index,load in enumerate(lbd):
            un.append(np.max(np.abs(disp[index])))

        # plot load-displacement curve - TODO: move to output / post-processing
        cnt = plot(un, lbd)

    out = min(step+1, abs(int(out_disp)))
    if out_disp > 0:
        u_out = un[out]
        l_out = lbd[out]
        print("\n************************************************************\n")
        print("Step: {0:2d} Load level: {1:.3f} Displacement: {2:.4e}".format(out, l_out,
                                                      u_out))
        print("\n************************************************************\n")
        return disp[out], sig, trac

    else:
        u_out = un[out]-un[out-1]
        l_out = lbd[out]-lbd[out-1]
        print("\n************************************************************\n")
        print("Step: {0:2d} Load level increment: {1:.3f} Displacement increment: {2:.4e}".format(out, l_out,
                                                      u_out))
        print("\n************************************************************\n")
        return disp[out]-disp[out-1], sig, trac


# plot the load-deflection curve
def plot(un, lbd):
    class Index(object):
        def stop(self, event):
            self.cnt = False
            plt.close()
        def add(self, event):
            self.cnt = True
            plt.close()

    callback = Index()
    callback.cnt=False
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.2)
    ax.plot(un, lbd, '-ok', color='black')
    ax.set(xlabel='displacement [mm]', ylabel='load factor [-]',
           title='')
    ax.grid()
    axstop = plt.axes([0.7, 0.05, 0.1, 0.075])
    axadd = plt.axes([0.81, 0.05, 0.1, 0.075])
    bstop = Button(axstop, 'stop')
    bstop.on_clicked(callback.stop)
    badd = Button(axadd, 'add')
    badd.on_clicked(callback.add)
    # fig.savefig("test.png")
    plt.show()

    return callback.cnt

# modify the global stiffness matrix and load vector for displacement boundary conditions
def bcGSM(gsm, glv, dispfaces):
    dim=len(glv)
    zero=np.zeros((dim), dtype=np.float64)
    dis=np.asarray(dispfaces)
    # fixdof=1: DOF is free; fixdof=0: DOF is fixed - used in calculation of residual load
    fixdof= np.ones((dim), dtype=int)

    for lf in dis:
        lx = lf[1][0] # True: free x-DOF; False: fixed x-DOF
        ly = lf[1][1] # True: free y-DOF; False: fixed y-DOF
        lz = lf[1][2] # True: free z-DOF; False: fixed z-DOF
        ux = uy = uz = 0.0

        if not lx : ux = lf[2][0] # prescribed displacement in x-direction
        if not ly : uy = lf[2][1] # prescribed displacement in y-direction
        if not lz : uz = lf[2][2] # prescribed displacement in z-direction

        for node in lf[0]:
            n3 = 3*(int(node)-1)
            if not lx:
                fixdof[n3]=0
                glv-=ux*gsm[n3]
                gsm[:,n3] = zero
                gsm[n3] = zero
                gsm[n3,n3] = 1.0
                glv[n3] = ux
            if not ly:
                fixdof[n3+1]=0
                glv-=uy*gsm[n3+1]
                gsm[:,n3+1] = zero
                gsm[n3+1] = zero
                gsm[n3+1,n3+1] = 1.0
                glv[n3+1] = uy
            if not lz:
                fixdof[n3+2]=0
                glv-=uz*gsm[n3+2]
                gsm[:,n3+2] = zero
                gsm[n3+2] = zero
                gsm[n3+2,n3+2] = 1.0
                glv[n3+2] = uz

    return gsm, glv, fixdof

# integrate stress-strain relationship for displacement increment du
def update_stress(elnodes, nocoord, materialbyElement, sig_yield, sig, du):

    u10 = np.zeros((30), dtype=np.float64) # displacements for the 10 tetrahedral nodes
    gp10, gp6 = gaussPoints()
    sig_update = np.array(np.zeros(24*len(elnodes), dtype=np.float64))
    pvec= np.array([1.0,1.0,1.0,0.0,0.0,0.0])


    for el, nodes in enumerate(elnodes):
        elpos=24*el
        dmat = hooke(el+1, materialbyElement)
        for index, nd in enumerate(nodes):
            n3=3*(nd-1)
            i3=3*index
            u10[i3]   = du[n3]
            u10[i3+1] = du[n3+1]
            u10[i3+2] = du[n3+2]
        xl = np.array([nocoord[nd-1] for nd in nodes]).T
        for index, ip in enumerate(gp10):
            ippos=elpos+6*index
            xi = ip[0]
            et = ip[1]
            ze = ip[2]
            xsj, shp, dshp, bmat = shape10tet(xi, et, ze, xl)

            # elastic test stress
            sig_test=sig[ippos:ippos+6]+np.dot(dmat,np.dot(bmat,u10))
            # von Mises stress
            normal=sig_test[:3]
            shear=sig_test[3:]
            pressure=np.average(normal)
            sig_mises=np.sqrt(1.5*np.linalg.norm(normal-pressure)**2+3.0*np.linalg.norm(shear)**2)
            # radial stress return to yield surface
            fac=np.minimum(sig_yield/sig_mises,1.0)
            sig_update[ippos:ippos+6] = fac * np.append(normal-pressure,shear) + pressure*pvec

    return np.array([sig_update])


def update_traction(interface_elements, nocoord, trac, du, kmax, kn, ks, shr_yield):
    u12 = np.zeros((36), dtype=np.float64) # displacements for the 10 tetrahedral nodes
    gp10, gp6 = gaussPoints()
    np10, np6 = nodalPoints()
    trac_update = np.array(np.zeros(18*len(interface_elements), dtype=np.float64))

    for el, nodes in enumerate(interface_elements):
        elpos = 18*el
        # material matrix in local coordinate system
        dmatloc = np.diag([kn*kmax,ks*kmax,ks*kmax])
        # coordinates of element nodes
        xl = np.array([nocoord[nd-1] for nd in nodes[:6]]).T
        # element nodal displacements
        for index, nd in enumerate(nodes):
            n3=3*(nd-1)
            i3=3*index
            u12[i3]   = du[n3]
            u12[i3+1] = du[n3+1]
            u12[i3+2] = du[n3+2]
        # contact stresses in the nodes (np6: Newton Cotes, gp6: Gauss)
        for index, ip in enumerate(np6):
            nppos = elpos + 3 * index
            xi = ip[0]
            et = ip[1]
            xsj, shp, bmat, xx, xt, xp = shape6tri(xi, et, xl)
            T=np.array([xp, xx, xt])
            dmatglob=np.dot(T.T,np.dot(dmatloc,T))
            # elastic test stress in local coordinates
            trac_test=trac[nppos:nppos+3] + np.dot(dmatglob,np.dot(bmat,u12))
            # contact stresses in local coordinates
            # TODO: bring stresses back to yield surface
            npstress = trac_test
            # add local contact stresses to global vectors
            trac_update[nppos:nppos+3]=npstress

    return np.array([trac_update])


# update internal load vector
def update_load(elnodes, interface_elements, nocoord, sig, trac):

    gp10, gp6 = gaussPoints()
    np10, np6 = nodalPoints()
    qin=np.array(np.zeros(3*len(nocoord), dtype=np.float64)) # internal load vector

    # volume element contribution
    for el, nodes in enumerate(elnodes):
        elv = np.zeros((30), dtype=np.float64)
        elpos = 24 * el
        xl = np.array([nocoord[nd-1] for nd in nodes]).T
        for index, ip in enumerate(gp10):
            xi = ip[0]
            et = ip[1]
            ze = ip[2]
            xsj, shp, dshp, bmat = shape10tet(xi, et, ze, xl)
            ippos=elpos+6*index
            elv += np.dot(bmat.T, sig[ippos:ippos+6]) * ip[3] * abs(xsj)
        for i in range(10):
            iglob = nodes[i] - 1
            iglob3 = 3 * iglob
            i3 = 3 * i
            for k in range(3):
                qin[iglob3 + k] += elv[i3 + k]

    # interface element contribution (np6: Newton Cotes, gp6: Gauss)
    for inel, nodes in enumerate(interface_elements):
        inelv = np.zeros((36), dtype=np.float64)
        elpos = 18 * inel
        xl = np.array([nocoord[nd-1] for nd in nodes[:6]]).T
        for index, ip in enumerate(np6):
            nppos = elpos + 3 * index
            xi = ip[0]
            et = ip[1]
            xsj, shp, bmat, xx, xt, xp = shape6tri(xi, et, xl)
            inelv += np.dot(bmat.T, trac[nppos:nppos+3]) * abs(xsj) * ip[2]
        for i in range(12):
            iglob = nodes[i] - 1
            iglob3 = 3 * iglob
            i3 = 3 * i
            for k in range (3):
                qin[iglob3 + k] += inelv[i3 + k]

    return qin


# map stresses to nodes
def mapStresses(elnodes, nocoord, interface_elements, disp, sig, trac, noce):

    # map maps corner node stresses to all tet10 nodes
    map = np.array([[1.0, 0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0],
                    [0.5, 0.5, 0.0, 0.0],
                    [0.0, 0.5, 0.5, 0.0],
                    [0.5, 0.0, 0.5, 0.0],
                    [0.5, 0.0, 0.0, 0.5],
                    [0.0, 0.5, 0.0, 0.5],
                    [0.0, 0.0, 0.5, 0.5]])

    expm = np.zeros((4, 4), dtype=np.float64) # extroplation matrix from Gauss points to corner nodes
    expm_int = np.zeros((6, 6), dtype=np.float64) # extroplation matrix from Integration Points to 6 tri6 nodes
    ipstress = np.zeros((4, 6), dtype=np.float64) # Tet10 stresses by Gauss point
    iptrac = np.zeros((6, 3), dtype=np.float64) # Tri6 tractions by integration point

    ip10, ip6 = gaussPoints()
    np10, np6 = nodalPoints()

    tet10stress=np.zeros((len(nocoord),6), dtype=np.float64)
    contactpressurevector=np.zeros((len(nocoord),3), dtype=np.float64)
    contactpressurevalue=np.zeros((len(nocoord)), dtype=np.float64)
    contactshearvector=np.zeros((len(nocoord),3), dtype=np.float64)

    xp_node=np.zeros((6, 3), dtype=np.float64) # normal vector in each of the 6 integration points
    xx_node=np.zeros((6, 3), dtype=np.float64) # Xi tangential vector in each of the 6 integration points
    xt_node=np.zeros((6, 3), dtype=np.float64) # shear vector ppd to the above 2 vectors in each of the 6 integration points

    step = len(sig)-1 # last step in the results

    # map stresses in volumes to nodal points
    for el, nodes in enumerate(elnodes):
        elpos=24*el
        xl = np.array([nocoord[nd-1] for nd in nodes]).T
        for index, ip in enumerate(ip10):
            xi = ip[0]
            et = ip[1]
            ze = ip[2]
            shp = shape4tet(xi, et, ze, xl)
            ippos=elpos+6*index
            ipstress[index]=sig[step][ippos:ippos+6] # ipstress (4x6): 6 stress components for 4 integration points
            for i in range (4):
                expm[index,i]=shp[i]
        expm_inv=np.linalg.inv(expm)
        npstress4 = np.dot(expm_inv,ipstress) # npstress4 (4x6): for each corner node (4) all stress components (6)
        numnodes = np.array([noce[nodes[n]-1] for n in range(10)]) # numnodes = number of nodes connected to node "nodes[n]-1"
        npstress10 = np.divide(np.dot(map, npstress4).T, numnodes).T # nodal point stress all nodes divided by number of connecting elements
        for index, nd in enumerate(nodes): tet10stress[nd-1] += npstress10[index]

    # For each interface element map the tractions to element nodes (np6: Newton Cotes, gp6: Gauss)
    # TODO: for extrapolated Gauss point results nodal averaging is required
    for el, nodes in enumerate(interface_elements):
        elpos=18*el
        xl = np.array([nocoord[nd-1] for nd in nodes[:6]]).T
        for index, ip in enumerate(np6):
            xi = ip[0]
            et = ip[1]
            xsj, shp, bmat, xx, xt, xp = shape6tri(xi, et, xl)
            xp_node[index] = xp
            xx_node[index] = xx
            xt_node[index] = xt
            T=np.array([xp, xx, xt])
            ippos = elpos + 3 * index
            iptrac[index] = np.dot(T, trac[step][ippos:ippos+3]) # local stresses at integration points
            for i in range (6):
                expm_int[index,i]=shp[i]
        expm_int_inv=np.linalg.inv(expm_int)
        nptrac = np.dot(expm_int_inv, iptrac) # local stresses extrapolated to nodes
        # add local contact stresses to global vectors
        for index, nd in enumerate(nodes[:6]):
            contactpressurevector[nd-1]=nptrac[index][0]*xp_node[index]
            contactpressurevalue[nd-1]=nptrac[index][0]
            contactshearvector[nd-1]=nptrac[index][1]*xx_node[index]+nptrac[index][2]*xt_node[index]

    return tet10stress, contactpressurevector, contactpressurevalue, contactshearvector

# fill resultobject with results
def pasteResults(doc, elnodes, nocoord, interface_elements, dis, tet10stress, contactpressurevector, contactpressurevalue, contactshearvector):

    analysis = doc.getObject("Analysis")

    if analysis == None:
        print("No Analysis object. Please create one first")
        raise SystemExit ()

    resVol = analysis.addObject(ObjectsFem.makeResultMechanical(doc))[0]
    resInt = analysis.addObject(ObjectsFem.makeResultMechanical(doc))[0]

    # VOLUME MESH START
    volnodes={}
    mode_disp_vol = {}
    elements_tetra10={}
    mode_results_vol = {}
    results=[]

    for index, coord in enumerate(nocoord):
        n3=3*index
        volnodes[index+1] = App.Vector(coord[0], coord[1], coord[2])
        mode_disp_vol[index+1] = App.Vector(dis[n3], dis[n3+1], dis[n3+2])

    for index, elem in enumerate(elnodes):
        elements_tetra10[index+1] = (elem[0],elem[2],elem[1],elem[3],elem[6],elem[5],elem[4],elem[7],elem[9],elem[8])

    mode_results_vol['disp'] = mode_disp_vol

    results.append(mode_results_vol)

    mvol = {
        'Nodes': volnodes,
        'Seg2Elem': {},
        'Seg3Elem': {},
        'Tria3Elem': {},
        'Tria6Elem': {},
        'Quad4Elem': {},
        'Quad8Elem': {},
        'Tetra4Elem': {},
        'Tetra10Elem': elements_tetra10,
        'Hexa8Elem': {},
        'Hexa20Elem': {},
        'Penta6Elem': {},
        'Penta15Elem': {},
        'Results': results
    }

    meshvol = itf.make_femmesh(mvol)

    result_mesh_object_1 = ObjectsFem.makeMeshResult(doc,'Result_Mesh_Volume')
    result_mesh_object_1.FemMesh = meshvol

    numnodes = len(nocoord)
    resVol.DisplacementVectors = [App.Vector(dis[3*n], dis[3*n+1], dis[3*n+2]) for n in range(numnodes)]
    resVol.DisplacementLengths = [np.linalg.norm([dis[3*n], dis[3*n+1], dis[3*n+2]]) for n in range(numnodes)]
    resVol.NodeStressXX = tet10stress.T[0].T.tolist()
    resVol.NodeStressYY = tet10stress.T[1].T.tolist()
    resVol.NodeStressZZ = tet10stress.T[2].T.tolist()
    resVol.NodeStressXY = tet10stress.T[3].T.tolist()
    resVol.NodeStressXZ = tet10stress.T[4].T.tolist()
    resVol.NodeStressYZ = tet10stress.T[5].T.tolist()

    resVol.Mesh = result_mesh_object_1
    resVol.NodeNumbers = [int(key) for key in resVol.Mesh.FemMesh.Nodes.keys()]

    resVol = itf.fill_femresult_mechanical(resVol, results)

    # VOLUME MESH FINISH

    # INTERFACE MESH START
    if interface_elements != []:

        intnodes={}
        mode_disp_int = {}
        intconnect={}
        newnode={}
        oldnode={}
        elements_tria6={}
        mode_results_int = {}
        results=[]

        index=0
        for i, intel in enumerate(interface_elements):
            for nd in intel[:6]:
                if nd not in intconnect:
                    index+=1
                    intconnect[nd]=index
                    intnodes[index] = App.Vector(nocoord[nd-1][0], nocoord[nd-1][1], nocoord[nd-1][2])
                    mode_disp_int[index] = App.Vector(contactpressurevector[nd-1][0], contactpressurevector[nd-1][1], contactpressurevector[nd-1][2])
                    newnode[nd]=index
                    oldnode[index]=nd

            elements_tria6[i + 1] = (newnode[intel[0]], newnode[intel[1]], newnode[intel[2]], newnode[intel[3]],
                                        newnode[intel[4]], newnode[intel[5]])


        mode_results_int['disp'] = mode_disp_int

        results.append(mode_results_int)

        mint = {
            'Nodes': intnodes,
            'Seg2Elem': {},
            'Seg3Elem': {},
            'Tria3Elem': {},
            'Tria6Elem': elements_tria6,
            'Quad4Elem': {},
            'Quad8Elem': {},
            'Tetra4Elem': {},
            'Tetra10Elem': {},
            'Hexa8Elem': {},
            'Hexa20Elem': {},
            'Penta6Elem': {},
            'Penta15Elem': {},
            'Results': results
        }

        meshint = itf.make_femmesh(mint)

        result_mesh_object_2 = ObjectsFem.makeMeshResult(doc,'Result_Mesh_Interface')
        result_mesh_object_2.FemMesh = meshint

        resInt.DisplacementVectors = [App.Vector(dis[3*(oldnode[n+1]-1)], dis[3*(oldnode[n+1]-1)+1], dis[3*(oldnode[n+1]-1)+2]) for n in range(len(oldnode))]
        resInt.DisplacementLengths = [contactpressurevalue[nd-1] for nd in oldnode.values()] # TODO: This is a dirty hack. move contact pressure to its own result object attribute

        resInt.Mesh = result_mesh_object_2
        resInt.NodeNumbers = [int(key) for key in resInt.Mesh.FemMesh.Nodes.keys()]

        resInt = itf.fill_femresult_mechanical(resInt, results)

    # INTERFACE MESH FINISH

    doc.recompute()

    return resInt, resVol