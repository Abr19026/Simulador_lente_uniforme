from collections import namedtuple
from dataclasses import dataclass
from typing import List
from itertools import chain

from math import asin, inf, sin, pi, isclose
from Classes import Ray, Collision, Vector, Point, Triangle
from lensDesigner import design_fresnel_prism_half, make_full_fresnel_lens


@dataclass(slots=True)
class Medium:
    name: str
    n: float

@dataclass
class RayMedium(Ray):
    medium: Medium

@dataclass
class Prism:
    triangle: Triangle
    medium: Medium
    def collision_with_ray(self, ray) -> List[Collision]:
        lens_inters = self.triangle.RayIntersections(ray)
        intersect = min(lens_inters, key=lambda inters: inters.distance, default=None)
        if intersect is None:
            return []
        normal = intersect.normal if intersect.from_outside else intersect.normal * -1
        return [Collision(normal, intersect.point, intersect.distance, intersect.from_outside, self, len(lens_inters) > 1)]

    def __repr__(self) -> str:
        return f"{self.triangle}"
        
@dataclass
class RayDistance(RayMedium):
    distance: float
    def from_medium(ray: RayMedium, distance):
        return RayDistance(ray.origin, ray.direction, ray.medium, distance)

    def __repr__(self) -> str:
        return f"origin: {self.origin.__repr__()}, dir: {self.direction.__repr__()} distance: {self.distance:,.3f}, end: {self.origin + self.direction.comps * self.distance}\n"
    
def find_min_collisions(ray, lenses: List[Prism]) -> List[Collision]:
    collisions = list(chain.from_iterable(map(lambda lens: lens.collision_with_ray(ray), lenses)))
    min_collision = min(collisions, key = lambda x: x.distance, default=None)
    if min_collision == None:
        return []
    return list(filter(lambda x: isclose(x.distance, min_collision.distance), collisions))


def refraction_vector(ray: Ray, collision_normal: Vector, ray_medium_n: float, new_medium_n: float) -> Vector:
    ray_comps = ray.get_unit_vec().comps
    # Very messy math, there are probably much cleaner and efficient ways to do it
    # First we get the angle from normal to ray by arcsin(y component of the ray vector on coordinates where normal vector is (1,0))
    thetha_1 = asin(collision_normal.comps.x * ray_comps.y - collision_normal.comps.y * ray_comps.x)
    try:
        # We apply snell's law
        thetha_2 = asin((ray_medium_n/new_medium_n) * sin(thetha_1))
        return ray.get_unit_vec().rotate(thetha_1 - thetha_2)
    except ValueError:
        # In case there is Total internal reflection (when argument to arcsin has abs. value greater than 1)
        return None


def next_ray(ray: RayMedium, lenses: List[Prism], outside_medium: Medium) -> RayMedium | None:
    # Obtiene colisiónes mas cercanas
    collisions = find_min_collisions(ray, lenses)
    # Finaliza rayo
    if len(collisions) < 1:
        return None
    
    collision = None
    if ray.medium == outside_medium:
        collision = next(filter(lambda x: x.from_outside and x.atraviesa, collisions))
    else:
        collision = next(filter(lambda x: not x.from_outside, collisions))
    
    new_medium = collision.col_object.medium if collision.from_outside else outside_medium
    # Crea rayo refractado
    dir_refraccion = refraction_vector(ray, collision.normal, ray.medium.n, new_medium.n)
    if dir_refraccion is None:
       return None
    return RayMedium(collision.pos, dir_refraccion, new_medium)

    

# Outside medium is air
def ray_trace(initial_rays: List[RayMedium], lenses: List[Prism], outside_medium: Medium) -> List[List[RayDistance]]:
    final_rays = []
    initial_rays = initial_rays.copy()
    # Por cada rayo
    for ray in initial_rays:
        pending_ray = ray
        ray_path = []
        while True:
            new_ray = next_ray(pending_ray, lenses, outside_medium)
            if new_ray is None:
                ray_path.append(RayDistance.from_medium(pending_ray, inf))
                break
            ray_path.append(RayDistance.from_medium(pending_ray, pending_ray.origin.distance(new_ray.origin)))
            pending_ray = new_ray
        final_rays.append(ray_path)
    return final_rays


def get_prism_data(theta, medium: Medium, a, b, half_N, r):
    half_lens = design_fresnel_prism_half(theta,medium.n,a,b,half_N,r)
    full_lens = make_full_fresnel_lens(half_lens)
    prism_data = []
    for prism in full_lens:
        triangle_points = tuple(Point(x[0], x[1]) for x in prism[0:3])
        prism_data.append(Prism(Triangle(triangle_points), medium))
    return prism_data


# NOTE:
# Program does not support rays shorter than epsilon units (Constant in Classes.py)
# Lenses must not overlap each other
if __name__ == "__main__":
    # Sets up the simulation environment
    medium_data = {"air": Medium("air", 1), "glass": Medium("glass", 1.5)}
    initial_rays = [RayMedium(Point(0, 0), Vector.init_from_angle(-13 * pi/180), medium_data["air"])]
    lenses = get_prism_data(25 * pi / 180, medium_data["glass"], 1, 4, 10, 4)
    final_rays = ray_trace(initial_rays,lenses, medium_data["air"])
    print(f"RAYS\n====\n{final_rays}")
    #print(f"\nLENSES\n======\n{lenses}")