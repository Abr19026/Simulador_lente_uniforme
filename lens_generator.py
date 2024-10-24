from math import sin, tan, cos, asin, atan, sqrt, pi

# theta: Medio ángulo del haz de luz
# n: índice de refracción de los prismas
# a: distancia horizontal desde la fuente al lente
# b: distancia horizontal del lente a la fuente
# N: Mitad de la cantidad de prismas del lente
# r: radio del disco
# Retorna una lista de tuplas-4 por prisma, cada tupla-4 contiene:
#   3 tuplas-2 que representan las coordenadas de los puntos en el prisma
#   El ángulo de apertura del prisma en radianes
def design_fresnel_prism_half(theta,n,a,b,N,r):
    # Other parameters
    max_beta_angle = 2 * pi / 5
    # Start 
    prism_data = []
    w = theta/N
    # for each prism
    for i in range(1, N+1):
        phi_0 = w * i               # Entrance angle
        phi_1 = asin(sin(phi_0)/n)  # first refraction angle
        h = a * tan(phi_0)          # y coordinate of the top of prism
        h_low = a * tan(phi_0 - w)  # y coordinate of bottom of prism
        f_height = r * sin(phi_0) / sin(theta) # height of the focus point
        th_out = atan((h - f_height)/b)
        fbeta = lambda beta: f_bnewton(beta, th_out,n,phi_1)
        dfbeta = lambda beta: df_bnewton(beta,th_out,n,phi_1)
        beta = Newton_method(fbeta, dfbeta, pi/4, 0.0005) # prism aperture
        if abs(beta) > max_beta_angle or abs(beta) + abs(th_out) > pi/2:
            print(f"angulo sobrepasado beta = {beta * 180 / pi}°, th_out = {th_out * 180 / pi}°")
            break
        prism_width = abs(tan(beta) * (h - h_low))
        prism_data.append(( (a,h_low), (a,h), (a + prism_width, h_low if beta > 0 else h), beta))
    return prism_data

# Takes as input the half lens created with design_fresnel_prism_half
def make_full_fresnel_lens(half_lens):
    full_lens = []
    for lens in half_lens:
        # Creates a vertically inverted copy of the prism
        inverted_copy = ((lens[0][0], -lens[0][1]), (lens[1][0], -lens[1][1]), (lens[2][0], -lens[2][1]), -lens[3])
        # Adds both the original lens and its inverted copy
        full_lens.append(lens)
        full_lens.append(inverted_copy)
    return full_lens

def Newton_method(f, df, initial_guess, error):
    x = initial_guess
    while abs(f(x)) > error:
        x = x - f(x) / df(x)
    return x

# function that becomes zero when the correct input (beta = lens aperture) is given
# t_out, n, phi_1 are parameters
def f_bnewton(beta, t_out, n, phi_1):
    return asin(sin(beta + t_out)/n) + phi_1 - beta

# derivative of f_b
def df_bnewton(beta, t_out, n, phi_1):
    return cos(beta + t_out) / (n * sqrt(1 - (sin(beta + t_out) ** 2) / (n ** 2))) - 1


#if __name__ == "__main__":
#    prismas = design_fresnel_prism_half(25 * pi / 180, 1.5, 1, 4, 10, 4)
#    for prisma in prismas:
#        print(prisma[3] * 180 / pi)
#        angle = atan((prisma[0][0] - prisma[2][0]) / (prisma[0][1] - prisma[1][1]))
#        if not isclose(angle * (-1 if prisma[3] < 0 else 1), prisma[3]):
#            print("No coinciden")