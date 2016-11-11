#!/usr/bin/env python

import math, sys
from PyQt4 import QtGui, QtCore

def makeVector(p1, p2, norm=False):
    return Vector(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z, norm)

def loadObj(filename):
    verts = []
    polys = []

    # read in each line, sort it appropriately as a vertex or face
    for line in open(filename, "r"):
        vals = line.split()

        if vals[0] == "v":
            v = Vector(float(vals[3]), float(vals[1]), float(vals[2]))
            verts.append(v)

        if vals[0] == "f":
            p1 = verts[int(vals[1])-1]
            p2 = verts[int(vals[2])-1]
            p3 = verts[int(vals[3])-1]

            p = Polygon(p1, p2, p3)
            polys.append(p)

    return polys

class Render(QtGui.QWidget):
    def __init__(self, scene, cam):
        super(Render, self).__init__()

        self.scene = scene
        self.cam = cam

        self.setGeometry(100, 100, self.cam.resX, self.cam.resY)
        self.setWindowTitle('Render')
        self.show()

    def paintEvent(self, e):
        qp = QtGui.QPainter()
        qp.begin(self)

        for y in range(self.cam.resY):
            for x in range(self.cam.resX):          
                # determine a ray through the current screen pixel
                ray = self.cam.pixel2ray((x, y))

                # trace that ray and return a luminance
                color = self.scene.trace(ray)

                print 'Solved', x, y, color
                
                qp.setPen(QtGui.QColor(color, color, color, 255))
                qp.drawPoint(x, y)

        qp.end()

class Scene(object):
    def __init__(self, objs=None):
        self.objs = []

        self.light = Vector(-10, -10, 10)

        for obj in objs:
            self.add(obj)

    def __iter__(self):
        for obj in self.objs:
            yield obj

    def add(self, obj):
        obj.setID(len(self.objs))
        self.objs.append(obj)

    def trace(self, ray, d=0):
        # begin with default values
        pointColor = 0
        color = 0
        closestT = None
        closestObj = None

        # iterate through each object in the scene, finding the first intersection
        for obj in self:
            t = obj.intersection(ray)

            if t:
                if not closestT or t < closestT:
                    closestT = t
                    closestObj = obj

        # if the ray strikes anything, determine shading
        if closestObj:
            i = ray.pointAt(closestT)

            normalVector = closestObj.normalAt(i)

            lightVector = makeVector(i, self.light, True)

            lightCos = normalVector.dot(lightVector)

            if lightCos < 0:
                pointColor = 0
            else:
                pointColor = lightCos * 255

            color = pointColor

        return color

class Vector(object):
    def __init__(self, x, y, z, norm=False):
        l = math.sqrt(math.pow(x, 2) + math.pow(y, 2) + math.pow(z, 2))

        self.x = x / l if norm else x
        self.y = y / l if norm else y
        self.z = z / l if norm else z

    def __add__(self, v):
        return Vector(self.x + v.x, self.y + v.y, self.z + v.z)

    def __sub__(self, v):
        return Vector(self.x - v.x, self.y - v.y, self.z - v.z)

    def __mul__(self, k):
        return Vector(self.x * k, self.y * k, self.z * k)

    def dot(self, v):
        return (self.x * v.x) + (self.y * v.y) + (self.z * v.z)

    def cross(self, v):
        x = (self.y * v.z) - (self.z * v.y)
        y = (self.z * v.x) - (self.x * v.z)
        z = (self.x * v.y) - (self.y * v.x)

        return Vector(x, y, z)

    def norm(self):
        return Vector(self.x, self.y, self.z, True)

class Ray(object):
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def pointAt(self, t):
        return self.a + (self.b * t)

class SceneObject(object):
    def __init__(self):
        self.id = None

    def setID(self, id):
        self.id = id

class Polygon(SceneObject):
    def __init__(self, v1, v2, v3):
        SceneObject.__init__(self)
        
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3

    def intersection(self, ray):
        n = self.normalAt()

        # return None if ray is parallel to the polygon's plane
        den = n.dot(ray.b)

        if den == 0:
            return None

        # return None if intersection lies behind the camera
        ti = n.dot(self.v1 - ray.a) / den

        if ti < 0:
            return None

        # determine if the intersection lies within the polygon
        i = ray.pointAt(ti)

        u = self.v2 - self.v1
        v = self.v3 - self.v1
        w = i - self.v1

        # define dot products so they're only calculated once
        uu = u.dot(u)
        uv = u.dot(v)
        uw = u.dot(w)
        vv = v.dot(v)
        vw = v.dot(w)

        den = math.pow(uv, 2) - (uu * vv)

        # return None if the intersection is out of bounds in the S dimension
        s = ((uv * vw) - (vv * uw)) / den

        if s < 0 or s > 1:
            return None
        
        # return None if the intersection is out of bounds in the T or S+T dimensions
        t = ((uv * uw) - (uu * vw)) / den

        if t < 0 or s + t > 1:
            return None

        return ti

    def normalAt(self, p=None):
        return (self.v3 - self.v2).cross(self.v1 - self.v2).norm()

class Sphere(SceneObject):
    def __init__(self, c, r):
        SceneObject.__init__(self)

        self.x = c[0]
        self.y = c[1]
        self.z = c[2]
        self.r = r

    def intersection(self, ray):
        A = math.pow(ray.v.dx, 2) + math.pow(ray.v.dy, 2) + math.pow(ray.v.dz, 2)
        B = 2 * ((ray.v.dx * (ray.x - self.x)) + (ray.v.dy * (ray.y - self.y)) + (ray.v.dz * (ray.z - self.z)))
        C = math.pow((ray.x - self.x), 2) + math.pow((ray.y - self.y), 2) + math.pow((ray.z - self.z), 2) - math.pow(self.r, 2)

        discriminant = math.pow(B, 2) - (4 * A * C)

        if discriminant <= 0:
            return None
        
        t1 = ((-1 * B) + math.sqrt(discriminant)) / (2 * A)
        t2 = ((-1 * B) - math.sqrt(discriminant)) / (2 * A)
        
        if t1 > 0: ti = t1
        if t2 > 0 and t2 < t1: ti = t2

        n = self.normalAt(ti)

        return ti

    def normalAt(self, p):
        return Vector(p[0] - self.x, p[1] - self.y, p[2] - self.z, True)

class Camera(object):
    def __init__(self, pos, rot, ang, res):
        self.x = pos[0]
        self.y = pos[1]
        self.z = pos[2]

        self.pan = rot[0]
        self.tilt = rot[1]

        self.resX = res[0]
        self.resY = res[1]

        self.ang = ang

        # calculate distance to screen in pixel units
        self.dts = self.resX / (2 * math.tan(self.ang / 2))

    # takes screen (X, Y) coordinates, returns a Ray object through that pixel
    def pixel2ray(self, coords):
        # calculate world coordinates without pan, tilt, or dolly
        normX = self.dts
        normY = (self.resX / 2) - coords[0]
        normZ = (self.resY / 2) - coords[1]
        xyd = math.sqrt(math.pow(normX, 2) + math.pow(normY, 2))

        # convert to spherical coordinates, include pan and tilt
        azm = math.atan(normY / normX) + self.pan
        inc = math.atan(normZ / xyd) + self.tilt
        rad = math.sqrt(math.pow(xyd, 2) + math.pow(normZ, 2))

        # convert back to Cartesian coordinates, include dolly
        trueX = (rad * math.cos(azm) * math.cos(inc)) + self.x
        trueY = (rad * math.sin(azm) * math.cos(inc)) + self.y
        trueZ = (rad * math.sin(inc)) + self.z

        # formulate origin and destination points for the ray
        po = Vector(self.x, self.y, self.z)
        pd = Vector(trueX, trueY, trueZ)

        return Ray(po, pd - po)

pos = (-10, 0.01, 1)
rot = (math.radians(0), math.radians(0))
ang = math.radians(60)
res = (360, 360)

cam = Camera(pos, rot, ang, res)

polys = loadObj("soccerball.obj")

scene = Scene(polys)

app = QtGui.QApplication(sys.argv)
render = Render(scene, cam)
sys.exit(app.exec_())