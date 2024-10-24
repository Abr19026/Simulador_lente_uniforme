from dataclasses import dataclass
from typing import List
from itertools import chain

from math import asin, inf, sin, pi, isclose
from Geometry import Ray, Vector, Point, Triangle
from lens_generator import design_fresnel_prism_half, make_full_fresnel_lens

@dataclass
class Collision:
    normal: Vector
    pos: Point
    distance: float
    from_outside: bool
    col_object: any
    atraviesa: bool

@dataclass(slots=True)
class Medium:
    name: str
    n: float
    opaque: bool

@dataclass
class RayMedium(Ray):
    medium: Medium

epsilon = 0.00001

@dataclass
class Prism:
    triangle: Triangle
    medium: Medium
    def collision_with_ray(self, ray) -> List[Collision]:
        lens_inters = list(filter(lambda intr: intr.distance > epsilon, self.triangle.RayIntersections(ray)))
        intersect = min(lens_inters, key=lambda inters: inters.distance, default=None)
        if intersect is None:
            return []
        normal = intersect.normal if intersect.from_outside else intersect.normal * -1
        atraviesa = len(lens_inters) > 1 and abs(lens_inters[0].distance - lens_inters[1].distance) > epsilon
        return [Collision(normal, intersect.point, intersect.distance, intersect.from_outside, self, atraviesa)]

    def __repr__(self) -> str:
        return f"{self.triangle}"
        
@dataclass
class RayDistance(RayMedium):
    distance: float
    def from_medium(ray: RayMedium, distance):
        return RayDistance(ray.origin, ray.direction, ray.medium, distance)

    def __repr__(self) -> str:
        return f"origin: {self.origin.__repr__()}, dir: {self.direction.__repr__()} distance: {self.distance:,.3f}, end: {self.origin + self.direction.comps * self.distance}\n"



def find_min_collisions(ray, objects: List[Prism]) -> List[Collision]:
    collisions = list(chain.from_iterable(map(lambda opt_object: opt_object.collision_with_ray(ray), objects)))
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


# It either returns a new ray, no ray, or direction of last collision
def next_ray(ray: RayMedium, objects: List[Prism], outside_medium: Medium) -> RayMedium | Collision | None:
    # Obtiene colisiónes mas cercanas
    collisions = find_min_collisions(ray, objects)
    # Finaliza rayo
    if len(collisions) < 1:
        return True
    # Choca rayo si es opaco
    opacos = list(filter(lambda col: col.col_object.medium.opaque, collisions))
    if len(opacos) > 0:
        return opacos[0]
    
    collision = None
    try:
        if ray.medium == outside_medium:
            collision = next(filter(lambda x: x.from_outside and x.atraviesa, collisions))
        else:
            collision = next(filter(lambda x: not x.from_outside, collisions))
    except StopIteration:
        return RayMedium(collisions[0].pos, ray.direction, ray.medium)
    
    new_medium = collision.col_object.medium if collision.from_outside else outside_medium
    # Crea rayo refractado
    dir_refraccion = refraction_vector(ray, collision.normal, ray.medium.n, new_medium.n)
    if dir_refraccion is None:
       return collision
    return RayMedium(collision.pos, dir_refraccion, new_medium)

    

# Outside medium is air
def ray_trace(initial_rays: List[RayMedium], objects: List[Prism], outside_medium: Medium) -> List[List[RayDistance]]:
    final_rays = []
    initial_rays = initial_rays.copy()
    # Por cada rayo
    for ray in initial_rays:
        pending_ray = ray
        ray_path = []
        while True:
            ray_or_collision = next_ray(pending_ray, objects, outside_medium)
            if ray_or_collision is None:
                ray_path.append(RayDistance.from_medium(pending_ray, inf))
                break
            elif type(ray_or_collision) is Collision:
                ray_path.append(RayDistance.from_medium(pending_ray, ray_or_collision.distance))
                break
            ray_path.append(RayDistance.from_medium(pending_ray, pending_ray.origin.distance(ray_or_collision.origin)))
            pending_ray = ray_or_collision
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

# light_half_angle = angulo al que emite luz el lente (<90)
# prism_refr_index = indice de refracción de los prismas de lente
# a = distancia del foco al lente
# d = distancia del foco al disco (d > a)
# half_N = mitad de la cantidad de prismas en el lente
# r = radio del disco
# I = cantidad de rayos a simular (Tratar que no sea multiplo ni divisor de half_N)
def simulate(light_half_angle: float, prism_refr_index: float, a: float, d: float, half_N: int, r: float, I: int):
    # Sets up the simulation environment
    # Creates mediums
    medium_data = {
        "air": Medium("air", 1, False), 
        "glass": Medium("glass", prism_refr_index, False), 
        "opaque": Medium("opaque", 0, True)}
    
    # Creates rays
    initial_rays = [RayMedium(Point(0, 0), Vector.init_from_angle(light_half_angle * (1 - 2/I * i)), medium_data["air"]) for i in range(0,I)]
    # Creates lens
    objects = get_prism_data(light_half_angle, medium_data["glass"], a, d - a, half_N, r)
    # Creates disc
    objects.append(Prism(Triangle((Point(d,r), Point(d,-r), Point(d+1,0))), medium_data["opaque"]))
    # Simulates rays
    ray_paths = ray_trace(initial_rays, objects, medium_data["air"])
    # Returns results 
    return ray_paths, objects

# NOTE:
# Program does not support rays shorter than epsilon units (Constant in Classes.py)
# Lenses must not overlap each other
if __name__ == "__main__":
    final_rays, objects = simulate(25 * pi / 180, 1.5, 1, 5, 10, 4, 44)
    #print(f"\nLENSES\n======\n{objects}")
    # checar que cada rayo termine cerca de r * sin(angulo_rayo_inicial) / sin(light_half_angle)
    print(f"RAYS\n======")
    for path in final_rays:
        print(path)

