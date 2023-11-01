from hilbert import omega
from misc.dual import Polar


omega_for_dual = Polar.omega(vertices=[(1.4, -0.1),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (-0.2, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])


primal_omega_weird = omega.Omega(vertices=[(1.4, -0.1),
                                  (1.1, 0.8),
                                  (0.6, 1.1),
                                  (-0.2, 1.2),
                                  (-.5, 0.8),
                                  (-.8, 0.3),
                                  (-.7, -.2)])

primal_omega_simplex = omega.Omega()