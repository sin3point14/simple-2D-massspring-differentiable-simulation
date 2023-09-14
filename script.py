import numpy as np

#//////////////////////////
#--------------------------
#  (0,0) \     //////  
#         \    ------
#    l1,k1 \     /  (1,0.5)
#           \  /  l2,k2
#            O   m
#         (0.5, 1)  
l1 = 1
l2 = 1
k1 = 2
k2 = 2

g = 9.81

spring1_end = np.array([0, 0])
spring2_end = np.array([1, 0.5])

particle_mass = 1
particle_position = np.array([0.5, 1])

target_position = np.array([0.5, 2])

def spring_force(position, l, k, end):
    stretch = np.linalg.norm(position - end) - l
    unit = (position - end) / np.linalg.norm(position - end)
    force = -k * stretch * unit
    return force

def spring_force_jacobian(position, l, k, end):
    curr_length = np.linalg.norm(position - end)
    unit = (position - end) / curr_length
    J = -k * ( (1 - l / curr_length) * (np.eye(2) - np.outer(unit, unit)) + np.outer(unit, unit) )
    return J

# solve for equilibrium using newton's method
F = lambda x: spring_force(x, l1, k1, spring1_end) + spring_force(x, l2, k2, spring2_end) + np.array([0, particle_mass * g])
J = lambda x: spring_force_jacobian(x, l1, k1, spring1_end) + spring_force_jacobian(x, l2, k2, spring2_end)

# verify jacobians
# for i in range(10):
#     x = np.random.rand(2) * 10
#     deltaX = np.random.rand(2)
#     dx = deltaX * (10**-i)
#     f = F(x)
#     f2 = F(x + dx)
#     print(np.linalg.norm(((f2 - f) / (10**-i)) - J(x) @ deltaX))

def equilibrium(particle_position, F, J):
    while np.linalg.norm(F(particle_position)) > 1e-6:
        particle_position = particle_position + -np.linalg.solve(J(particle_position), F(particle_position))
    return particle_position

particle_position = equilibrium(particle_position, F, J)
print("Initial equilibrium:")
print("Positions: ", particle_position)


def delf_delk(position, l, end):
    curr_length = np.linalg.norm(position - end)
    unit = (position - end) / curr_length
    return -(curr_length - l) * unit

T = lambda x: np.linalg.norm(x - target_position)**2
# gradient descent
alpha = 0.1
while np.linalg.norm(T(particle_position)) > 1e-3:
    # delf_delx * delx_delp = - delf_delp
    # p are our control parameters: k1, k2 in this case
    delf_delp = np.column_stack((delf_delk(particle_position, l1, spring1_end), delf_delk(particle_position, l2, spring2_end)))
    delx_delp = -np.linalg.inv(J(particle_position)) @ delf_delp

    grad_T = (particle_position - target_position) @ delx_delp

    delta_p = - alpha * grad_T
    k1 += delta_p[0]
    k2 += delta_p[1]
    # recalculate equilibrium
    particle_position = equilibrium(particle_position, F, J)

print("Final equilibrium:")
print("Spring constants: ", k1, k2)
print("Positions: ", particle_position)

