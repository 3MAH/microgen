import os
import math
from numpy import cos, sin, pi, abs
import pygalmesh


class Custom(pygalmesh.DomainBase):
    def __init__(self, rve, function, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1
        self.function = function

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return eval(self.function) + self.height


class Hyperboloid(pygalmesh.DomainBase):
    def __init__(self, max_edge_size_at_feature_edges):
        super().__init__()
        self.z0 = -1.0
        self.z1 = 1.0
        self.waist_radius = 0.5
        self.max_edge_size_at_feature_edges = max_edge_size_at_feature_edges

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        if self.z0 < z and z < self.z1:
            return x ** 2 + y ** 2 - (z ** 2 + self.waist_radius) ** 2
        return 1.0

    def get_bounding_sphere_squared_radius(self):
        z_max = max(abs(self.z0), abs(self.z1))
        r_max = z_max**2 + self.waist_radius
        return r_max * r_max + z_max * z_max

    def getFeatures(self):
        radius0 = self.z0**2 + self.waist_radius
        n0 = int(2 * pi * radius0 / self.max_edge_size_at_feature_edges)
        circ0 = [
            [
                radius0 * cos((2 * pi * k) / n0),
                radius0 * sin((2 * pi * k) / n0),
                self.z0,
            ]
            for k in range(n0)
        ]
        circ0.append(circ0[0])

        radius1 = self.z1**2 + self.waist_radius
        n1 = int(2 * pi * radius1 / self.max_edge_size_at_feature_edges)
        circ1 = [
            [
                radius1 * cos((2 * pi * k) / n1),
                radius1 * sin((2 * pi * k) / n1),
                self.z1,
            ]
            for k in range(n1)
        ]
        circ1.append(circ1[0])
        return [circ0, circ1]


class SchwarzP(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z) + self.height


class SchwarzD(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        a = sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z)
        b = sin(2*pi*x) * cos(2*pi*y) * cos(2*pi*z)
        c = cos(2*pi*x) * sin(2*pi*y) * cos(2*pi*z)
        d = cos(2*pi*x) * cos(2*pi*y) * sin(2*pi*z)
        return a + b + c + d + self.height


class Neovius(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]

        a = 3*cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)
        b = 4*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)

        return a + b + self.height

class SchoenIWP(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]

        l = 2*(cos(2*pi*x)*cos(2*pi*y) + \
            cos(2*pi*y)*cos(2*pi*z) + \
            cos(2*pi*z)*cos(2*pi*x))
        m = cos(4*pi*x) + cos(4*pi*y) + cos(4*pi*z)
        
        return l - m + self.height
    
class SchoenFRD(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]

        a = 4*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)
        b = cos(4*pi*x)*cos(4*pi*y) + \
            cos(4*pi*y)*cos(4*pi*z) + \
            cos(4*pi*z)*cos(4*pi*x)
        return a - b + self.height
    
class FischerKochS(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]

        a = cos(4*pi*x)*sin(2*pi*y)*cos(2*pi*z)
        b = cos(2*pi*x)*cos(4*pi*y)*sin(2*pi*z)
        c = sin(2*pi*x)*cos(2*pi*y)*cos(4*pi*z)

        return a + b + c + self.height
    
class PMY(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]

        a = 2*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)
        b = sin(4*pi*x) * sin(2*pi*y)
        c = sin(2*pi*x) * sin(4*pi*z)
        d = sin(4*pi*y) * sin(2*pi*z)

        return a + b + c + d + self.height
    
class HoneyComb(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return sin(2*pi*x)*cos(2*pi*y) + sin(2*pi*y) + cos(2*pi*z) + self.height
        
class Gyroid(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        if abs(x) + abs(y) + abs(z) > 1.0e-8:
            return (sin(2*pi*x)*cos(2*pi*y)
                   + sin(2*pi*y)*cos(2*pi*z) 
                   + sin(2*pi*z)*cos(2*pi*x)
                   + self.height)
        else:
            return 1.0


class Wolfram_P(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return sin(2*pi*x) + sin(2*pi*y) + sin(2*pi*z) + self.height


class Wolfram_iWP(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return cos(2*pi*x)*cos(2*pi*y) + cos(2*pi*y)*cos(2*pi*z) \
               + cos(2*pi*x)*cos(2*pi*z) \
               - 3 * cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z) \
               + self.height


class Wolfram_P_W(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return 4*(cos(2*pi*x)*cos(2*pi*y) \
                  + cos(2*pi*y)*cos(2*pi*z) \
                  + cos(2*pi*x)*cos(2*pi*z)) - \
               3*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z) + self.height


class Wolfram_lindinoid(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return 0.5*(sin(2*2*pi*x)*cos(2*pi*y)*sin(2*pi*z) \
                    + sin(2*2*pi*y)*cos(2*pi*z)*sin(2*pi*x) \
                    + sin(2*2*pi*z)*cos(2*pi*x)*sin(2*pi*y)) \
                - 0.5*(cos(2*2*pi*x)*cos(2*2*pi*y) \
                       + cos(2*2*pi*y)*cos(2*2*pi*z) \
                       + cos(2*2*pi*z)*cos(2*2*pi*x)) + .15 + self.height


class Wolfram_DoubleGyroid(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return 2.75*(sin(2*2*pi*x)*sin(2*pi*z)*cos(2*pi*y) \
                     + sin(2*2*pi*y)*sin(2*pi*x)*cos(2*pi*z) \
                     + sin(2*2*pi*z)*sin(2*pi*y)*cos(2*pi*x)) \
               - 1*(cos(2*2*pi*x)*cos(2*2*pi*y) \
                  + cos(2*2*pi*y)*cos(2*2*pi*z) \
                  + cos(2*2*pi*z)*cos(2*2*pi*x)) + self.height


class Wolfram_Gprime(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return (5*(sin(2*2*pi*x)*sin(2*pi*z)*cos(2*pi*y) 
                  + sin(2*2*pi*y)*sin(2*pi*x)*cos(2*pi*z)  
                  + sin(2*2*pi*z)*sin(2*pi*y)*cos(2*pi*x))  
               + 1*(cos(2*2*pi*x)*cos(2*2*pi*y) 
                  + cos(2*2*pi*y)*cos(2*2*pi*z) 
                  + cos(2*2*pi*z)*cos(2*2*pi*x)) + self.height)



class Wolfram_DoubleDiamond(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return (1*(sin(2*2*pi*x)*sin(2*2*pi*y) 
                  + sin(2*2*pi*y)*sin(2*2*pi*z) 
                  + sin(2*2*pi*x)*sin(2*2*pi*z)) 
               + cos(2*2*pi*x)*cos(2*2*pi*y)*cos(2*2*pi*z) 
               + self.height)


class Wolfram_Dprime(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return (sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z) 
               + cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z) 
               - 1*(cos(2*2*pi*x)*cos(2*2*pi*y) 
                    + cos(2*2*pi*y)*cos(2*2*pi*z) 
                    + cos(2*2*pi*z)*cos(2*2*pi*x)) - 0.4 + self.height)



class Wolfram_DoubleP(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return (0.5*(cos(2*pi*x)*cos(2*pi*y) 
                    + cos(2*pi*y)*cos(2*pi*z) 
                    + cos(2*pi*z)*cos(2*pi*x))  
               + 0.2*(cos(2*2*pi*x) 
                    + cos(2*2*pi*y) 
                    + cos(2*2*pi*z)) + self.height)


class Wolfram_OCTO(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return (4*(cos(2*pi*x)*cos(2*pi*y) 
                  + cos(2*pi*y)*cos(2*pi*z) 
                  + cos(2*pi*z)*cos(2*pi*x)) 
               - 2.8*(cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)) 
               + 1*(cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)) 
               + 1.5 + self.height)



class Wolfram_PN(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return (0.6*(cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)) 
               + 0.4*(cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)) 
               + 0.2*(cos(2*2*pi*x)*cos(2*2*pi*y)*cos(2*2*pi*z)) 
               + 0.2*(cos(2*2*pi*x) + cos(2*2*pi*y) + cos(2*2*pi*z))   
               + 0.1*(cos(3*2*pi*x) + cos(3*2*pi*y) + cos(3*2*pi*z)) 
               + 0.2*(cos(2*pi*x)*cos(2*pi*y) 
                      + cos(2*pi*y)*cos(2*pi*z) 
                      + cos(2*pi*z)*cos(2*pi*x)) + self.height)


class Wolfram_KP(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return (0.6*(cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z))  
               + 0.7*(cos(2*pi*x)*cos(2*pi*y) 
                      + cos(2*pi*y)*cos(2*pi*z) 
                      + cos(2*pi*z)*cos(2*pi*x)) 
               - 0.9*(cos(2*2*pi*x)*cos(2*2*pi*y)*cos(2*2*pi*z)) 
               + 0.4 + self.height)


class Wolfram_FRD(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return (8*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z) 
               + 1*(cos(2*2*pi*x)*cos(2*2*pi*y)*cos(2*2*pi*z))  
               - 1*(cos(2*2*pi*x)*cos(2*2*pi*y) 
                    + cos(2*2*pi*y)*cos(2*2*pi*z) 
                    + cos(2*2*pi*z)*cos(2*2*pi*x)) 
               + self.height)



class Wolfram_splitP(pygalmesh.DomainBase):
    def __init__(self, rve, height):
        super().__init__()
        self.height = height
        self.z0 = 0.0
        self.z1 = rve.dz
        self.waist_radius = math.sqrt((0.5 * rve.dx) ** 2 + (0.5 * rve.dy) ** 2)
        self.bounding_sphere_squared_radius = math.sqrt((0.5 * rve.dx) ** 2 +
                                                        (0.5 * rve.dy) ** 2 +
                                                        (0.5 * rve.dz) ** 2) * 1.1

    def get_bounding_sphere_squared_radius(self):
        return self.bounding_sphere_squared_radius

    def eval(self, pos):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        return (1.1*(sin(2*2*pi*x)*sin(2*pi*z)*cos(2*pi*y) 
                    + sin(2*2*pi*y)*sin(2*pi*x)*cos(2*pi*z) 
                    + sin(2*2*pi*z)*sin(2*pi*y)*cos(2*pi*x)) 
               - 0.2*(cos(2*2*pi*x)*cos(2*2*pi*y) 
                      + cos(2*2*pi*y)*cos(2*2*pi*z) 
                      + cos(2*2*pi*z)*cos(2*2*pi*x)) 
               - 0.4*(cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)) 
               + self.height)


def generateTPMS(
    type_tpms,
    thickness,
    rve,
    sizeMesh=0.05,
    minFacetAngle=10.0,
    maxRadius=0.05,
    path_data="",
    function=""
):

    thickness = thickness * pi

    if type_tpms == "custom":
        s_testplus = Custom(rve, function, thickness / 4.0)
        s_testminus = Custom(rve, function, -thickness / 4.0)
        s_plus = Custom(rve, function, thickness / 2.0)
        s_minus = Custom(rve, function, -1.0 * thickness / 2.0)  
    elif type_tpms == "honeycomb":
        s_testplus = HoneyComb(rve, thickness / 4.0)
        s_testminus = HoneyComb(rve, -thickness / 4.0)
        s_plus = HoneyComb(rve, thickness / 2.0)
        s_minus = HoneyComb(rve, -1.0 * thickness / 2.0)    
    elif type_tpms == "gyroid":
        s_testplus = Gyroid(rve, thickness / 4.0)
        s_testminus = Gyroid(rve, -thickness / 4.0)
        s_plus = Gyroid(rve, thickness / 2.0)
        s_minus = Gyroid(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "schwarzP":
        s_testplus = SchwarzP(rve, thickness / 4.0)
        s_testminus = SchwarzP(rve, -thickness / 4.0)
        s_plus = SchwarzP(rve, thickness / 2.0)
        s_minus = SchwarzP(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "schwarzD":
        s_testplus = SchwarzD(rve, thickness / 4.0)
        s_testminus = SchwarzD(rve, -thickness / 4.0)
        s_plus = SchwarzD(rve, thickness / 2.0)
        s_minus = SchwarzD(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "neovius":
        s_testplus = Neovius(rve, thickness / 4.0)
        s_testminus = Neovius(rve, -thickness / 4.0)
        s_plus = Neovius(rve, thickness / 2.0)
        s_minus = Neovius(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "schoenIWP":
        s_testplus = SchoenIWP(rve, thickness / 4.0)
        s_testminus = SchoenIWP(rve, -thickness / 4.0)
        s_plus = SchoenIWP(rve, thickness / 2.0)
        s_minus = SchoenIWP(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "schoenFRD":
        s_testplus = SchoenFRD(rve, thickness / 4.0)
        s_testminus = SchoenFRD(rve, -thickness / 4.0)
        s_plus = SchoenFRD(rve, thickness / 2.0)
        s_minus = SchoenFRD(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "fischerKochS":
        s_testplus = FischerKochS(rve, thickness / 4.0)
        s_testminus = FischerKochS(rve, -thickness / 4.0)
        s_plus = FischerKochS(rve, thickness / 2.0)
        s_minus = FischerKochS(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "pmy":
        s_testplus = PMY(rve, thickness / 4.0)
        s_testminus = PMY(rve, -thickness / 4.0)
        s_plus = PMY(rve, thickness / 2.0)
        s_minus = PMY(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "Wolfram_P":
        s_testplus = Wolfram_P(rve, thickness / 4.0)
        s_testminus = Wolfram_P(rve, -thickness / 4.0)
        s_plus = Wolfram_P(rve, thickness / 2.0)
        s_minus = Wolfram_P(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "Wolfram_iWP":
        s_testplus = Wolfram_iWP(rve, thickness / 4.0)
        s_testminus = Wolfram_iWP(rve, -thickness / 4.0)
        s_plus = Wolfram_iWP(rve, thickness / 2.0)
        s_minus = Wolfram_iWP(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "Wolfram_P_W":
        s_testplus = Wolfram_P_W(rve, thickness / 4.0)
        s_testminus = Wolfram_P_W(rve, -thickness / 4.0)
        s_plus = Wolfram_P_W(rve, thickness / 2.0)
        s_minus = Wolfram_P_W(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "Wolfram_lindinoid":
        s_testplus = Wolfram_lindinoid(rve, thickness / 4.0)
        s_testminus = Wolfram_lindinoid(rve, -thickness / 4.0)
        s_plus = Wolfram_lindinoid(rve, thickness / 2.0)
        s_minus = Wolfram_lindinoid(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "Wolfram_DoubleGyroid":
        s_testplus = Wolfram_DoubleGyroid(rve, thickness / 4.0)
        s_testminus = Wolfram_DoubleGyroid(rve, -thickness / 4.0)
        s_plus = Wolfram_DoubleGyroid(rve, thickness / 2.0)
        s_minus = Wolfram_DoubleGyroid(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "Wolfram_Gprime":
        s_testplus = Wolfram_Gprime(rve, thickness / 4.0)
        s_testminus = Wolfram_Gprime(rve, -thickness / 4.0)
        s_plus = Wolfram_Gprime(rve, thickness / 2.0)
        s_minus = Wolfram_Gprime(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "Wolfram_DoubleDiamond":
        s_testplus = Wolfram_DoubleDiamond(rve, thickness / 4.0)
        s_testminus = Wolfram_DoubleDiamond(rve, -thickness / 4.0)
        s_plus = Wolfram_DoubleDiamond(rve, thickness / 2.0)
        s_minus = Wolfram_DoubleDiamond(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "Wolfram_Dprime":
        s_testplus = Wolfram_Dprime(rve, thickness / 4.0)
        s_testminus = Wolfram_Dprime(rve, -thickness / 4.0)
        s_plus = Wolfram_Dprime(rve, thickness / 2.0)
        s_minus = Wolfram_Dprime(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "Wolfram_DoubleP":
        s_testplus = Wolfram_DoubleP(rve, thickness / 4.0)
        s_testminus = Wolfram_DoubleP(rve, -thickness / 4.0)
        s_plus = Wolfram_DoubleP(rve, thickness / 2.0)
        s_minus = Wolfram_DoubleP(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "Wolfram_OCTO":
        s_testplus = Wolfram_OCTO(rve, thickness / 4.0)
        s_testminus = Wolfram_OCTO(rve, -thickness / 4.0)
        s_plus = Wolfram_OCTO(rve, thickness / 2.0)
        s_minus = Wolfram_OCTO(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "Wolfram_PN":
        s_testplus = Wolfram_PN(rve, thickness / 4.0)
        s_testminus = Wolfram_PN(rve, -thickness / 4.0)
        s_plus = Wolfram_PN(rve, thickness / 2.0)
        s_minus = Wolfram_PN(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "Wolfram_KP":
        s_testplus = Wolfram_KP(rve, thickness / 4.0)
        s_testminus = Wolfram_KP(rve, -thickness / 4.0)
        s_plus = Wolfram_KP(rve, thickness / 2.0)
        s_minus = Wolfram_KP(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "Wolfram_FRD":
        s_testplus = Wolfram_FRD(rve, thickness / 4.0)
        s_testminus = Wolfram_FRD(rve, -thickness / 4.0)
        s_plus = Wolfram_FRD(rve, thickness / 2.0)
        s_minus = Wolfram_FRD(rve, -1.0 * thickness / 2.0)
    elif type_tpms == "Wolfram_splitP":
        s_testplus = Wolfram_splitP(rve, thickness / 4.0)
        s_testminus = Wolfram_splitP(rve, -thickness / 4.0)
        s_plus = Wolfram_splitP(rve, thickness / 2.0)
        s_minus = Wolfram_splitP(rve, -1.0 * thickness / 2.0)
    else:
        print("Error, the tpms is not recognized")
        return False

    mesh_surf_testplus = pygalmesh.generate_surface_mesh(
        s_testplus,
        min_facet_angle=minFacetAngle,
        max_radius_surface_delaunay_ball=maxRadius,
        max_facet_distance=sizeMesh,
    )

    mesh_surf_testminus = pygalmesh.generate_surface_mesh(
        s_testminus,
        min_facet_angle=minFacetAngle,
        max_radius_surface_delaunay_ball=maxRadius,
        max_facet_distance=sizeMesh,
    )

    mesh_surf_plus = pygalmesh.generate_surface_mesh(
        s_plus,
        min_facet_angle=minFacetAngle,
        max_radius_surface_delaunay_ball=maxRadius,
        max_facet_distance=sizeMesh,
    )

    mesh_surf_minus = pygalmesh.generate_surface_mesh(
        s_minus,
        min_facet_angle=minFacetAngle,
        max_radius_surface_delaunay_ball=maxRadius,
        max_facet_distance=sizeMesh,
    )

    if path_data != '':
        if not os.path.isdir(path_data):
            os.mkdir(path_data)
        mesh_surf_testplus.write(path_data + '/' + 'tpms_testplus.stl')
        mesh_surf_testminus.write(path_data + '/' + 'tpms_testminus.stl')
        mesh_surf_plus.write(path_data + '/' + 'tpms_plus.stl')
        mesh_surf_minus.write(path_data + '/' + 'tpms_minus.stl')
    else:
        mesh_surf_testplus.write("tpms_testplus.stl")
        mesh_surf_testminus.write("tpms_testminus.stl")
        mesh_surf_plus.write("tpms_plus.stl")
        mesh_surf_minus.write("tpms_minus.stl")

    return True

def is_function_valid(function):
    try:
        x = 0
        y = 0
        z = 0
        eval(function)
    except NameError:
        return False
    except SyntaxError:
        return False
    return True
