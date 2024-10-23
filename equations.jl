#=

Hess-Smith Panel Method

Governing equation (Laplace) is automatically satisfied,
as is the farfield boundary condition.

Satisfy flow tangency and the Kutta condition

Place the source/vortex distributions on the surface of the body
rather than at the chord line, and use the exact 
flow tangency conditions imposed at the surface.

Model the flow = freestream, line source distributions
Eq.2.52 (source line distribution -- 
the potential function expressed as an integral along the line distribution)
and vortex line distributions Eq.2.59 
(line of vortices with strength per unit length -- potential function as an integral along this distribution).

For boundary conditions:
Choosing sources and vortices (singularities) automatically satisfies the condition that induced velocity must go to zero in the farfield.

We also need a no flow-through condition (total velocity vector is tangent to the airfoil surface (2.72 and 2.73).

v_t (2.84)
vortex line integrals: 2.85 and 2.86

We can know the source distribution simply from the known airfoil thickness distribution (2.89)

Vortex distribution is not as simple, and requires solving an integral equation with a known AoA and camber distribution. 2.90

=#