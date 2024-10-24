from collections import namedtuple
from typing import Tuple, List
from math import isclose, sqrt, sin, cos, acos, pi
from dataclasses import dataclass



@dataclass
class Point:
    x: float
    y: float
    
    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y)
    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y)
    def __mul__(self, factor):
        return Point(self.x * factor , self.y * factor)
    def __truediv__(self, factor: float):
        return Point(self.x / factor, self.y / factor)
    def midpoint(self, other):
        return (self + other) / 2
    def distance(self, other):
        return sqrt((self.x - other.x) ** 2 + (self.y - other.y) ** 2)
    def __repr__(self) -> str:
        return f"({self.x:,.4f}, {self.y:,.4f})"
        
@dataclass
class Vector:
    comps: Point
    def get_length(self):
        return self.comps.distance(Point(0,0))
    
    @staticmethod
    def init(x_comp, y_comp):
        return Vector(Point(x_comp, y_comp))

    def dot_product(self, other):
        return self.comps.x * other.comps.x  + self.comps.y * other.comps.y  
    
    @staticmethod
    def init_from_angle(angle):
        return Vector(Point(cos(angle), sin(angle)))

    def get_unit_vec(self):
        return Vector(self.comps / self.get_length())
    
    def rotate(self, rad_angle):
        return Vector(Point(cos(rad_angle) * self.comps.x - sin(rad_angle) * self.comps.y,
                            sin(rad_angle) * self.comps.x + cos(rad_angle) * self.comps.y))
    
    def get_abs_angle(self):
        angle = acos(self.get_unit_vec().comps.x)
        if self.comps.y < 0:
            angle = -angle
        return angle

    def __repr__(self):
        return f"[{self.comps.x:,.4f}, {self.comps.y:,.4f}] ({self.get_abs_angle() * 180 / pi:,.2f} deg)"
    
    def __mul__(self, factor):
        return Vector(self.comps * factor)


@dataclass
class Ray:
    origin: Point
    direction: Vector
    
    def get_unit_vec(self) -> Vector:
        return self.direction

@dataclass
class RayIntersection:
    point: Point
    normal: Vector
    distance: float
    from_outside: bool


@dataclass
class Triangle:
    points: Tuple[Point,Point,Point]

    def RayIntersections(self, ray: Ray) -> List[RayIntersection]:
        intersections = []
        ray_dir = ray.get_unit_vec()
        centro = sum(self.points, Point(0,0))/3
        for segment in ((self.points[0], self.points[1]),(self.points[0], self.points[2]),(self.points[1], self.points[2])):
            delta_x = segment[1].x - segment[0].x
            delta_y = segment[1].y - segment[0].y
            divisor = delta_y * ray_dir.comps.x - delta_x * ray_dir.comps.y
            if divisor == 0:
                continue
            distance_intersection = (delta_x * (ray.origin.y - segment[0].y) - delta_y * (ray.origin.x - segment[0].x)) / divisor
            inter_point: Point = ray.origin + ray_dir.comps * distance_intersection
            
            def in_float_range(value, min_v, max_v):
                return (value > min_v and value < max_v) or isclose(value,min_v) or isclose(value,max_v)

            #if distance_intersection > epsilon and \
            if in_float_range(inter_point.x, min(segment, key=lambda s:s.x).x, max(segment, key=lambda s:s.x).x) and \
            in_float_range(inter_point.y, min(segment, key=lambda s:s.y).y, max(segment, key=lambda s:s.y).y):
                normal = Vector(Point(delta_y, -delta_x))
                if normal.dot_product(Vector(centro-inter_point)) > 0:
                    normal = normal * -1
                from_outside = normal.dot_product(ray_dir) < 0
                intersections.append(RayIntersection(inter_point, normal.get_unit_vec(), distance_intersection, from_outside))
        
        return intersections

    def __repr__(self) -> str:
        return f"{self.points}"