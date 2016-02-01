import numpy as np

_InIceRelativeDepth = 1455 # IceCube depth wrt IceTop in m

#------------------------------------------------------------------------
#------------------------------------------------------------------------

#def containCore(core, vertices):
#  arrayEdge = mp.Path(np.array(vertices))
#  return arrayEdge.contains_point(core)

def containCore(core, poly):
    x = core[0]
    y = core[1]
    n = len(poly)
    inside =False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

#------------------------------------------------------------------------

def shrinkEdge(points, scale):
    sx = []
    sy = []

    for point in points:
      sx.append(point[0] * scale)
      sy.append(point[1] * scale)

    return (sx, sy)
      
#------------------------------------------------------------------------

def hitsInIce(core, dir, vertices, deltaZ = _InIceRelativeDepth):
  x = core[0]
  y = core[1]

  zen = dir[0]
  azi = dir[1]

  xi = - np.tan(zen) * np.cos(azi) * deltaZ + x
  yi = - np.tan(zen) * np.sin(azi) * deltaZ + y

  return containCore([xi, yi], vertices)

#------------------------------------------------------------------------

def distanceToInnerString(x, y, zen, azi, dom_x, dom_y, dom_z, deltaZ = _InIceRelativeDepth):

  xi = - np.tan(zen) * np.cos(azi) * (deltaZ-dom_z) + x
  yi = - np.tan(zen) * np.sin(azi) * (deltaZ-dom_z) + y

  distance = np.sqrt((xi-dom_x)**2 + (yi-dom_y)**2)
  return distance

#------------------------------------------------------------------------

def loadGeometryTxt(name):
    file = np.genfromtxt(name)
    x = file[:,0]
    y = file[:,1]
    return zip(x, y)

#------------------------------------------------------------------------

def coincidentCutDirPos(partDir, partPos, itgeo, icgeo):
  itdir = [partDir.dir.zenith, partDir.dir.azimuth]
  itpos = [partPos.pos.x, partPos.pos.y]
 
  ITCut = containCore(itpos, itgeo)
  ICCut = hitsInIce(itpos, itdir, icgeo)

  result = (ITCut & ICCut)
  return result
'''
def coincidentCutDirPos(itdir, itpos, itgeo, icgeo):
  ITCut = containCore(itpos, itgeo)
  ICCut = hitsInIce(itpos, itdir, icgeo)

  result = (ITCut & ICCut)
  return result
'''
#------------------------------------------------------------------------
def distanceToTrack(pulse_x, pulse_y, pulse_z, core, dir, deltaZ = _InIceRelativeDepth):
  x = core[0]
  y = core[1]

  zen = dir[0]
  azi = dir[1]

  xi = - np.tan(zen) * np.cos(azi) * (deltaZ + 500 - pulse_z) + x
  yi = - np.tan(zen) * np.sin(azi) * (deltaZ + 500 - pulse_z) + y

  d = np.sqrt((xi-pulse_x)**2 + (yi-pulse_y)**2)
  return d

#------------------------------------------------------------------------
