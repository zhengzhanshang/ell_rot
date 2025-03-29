#main code
# spin-up in a sphere
# start with: python spheroid.py input_spheroid.json

import json
import sys
import numpy as np
import dedalus.public as d3
import logging
import random
logger = logging.getLogger(__name__)


restart = (len(sys.argv) > 1 and sys.argv[1] == '--restart')

# Parameters
if len(sys.argv)==1:
   print('Argument fehlt!')
   exit()
else:
   params=json.load(open(sys.argv[1]))

Radius = params['Radius']
ell = params['ell']
a = params['a']
Nphi = params['Nphi']
Ntheta = params['Ntheta']
Nr = params['Nr']
Ek = params['Ek']

volume = 4*np.pi*Radius**3/(3*(1-ell)) #c=Radius

stop_sim_time = params['stop_sim_time'] + params['stop_sim_time']*restart
max_timestep = params['max_timestep']

run_name = params['run_name']
ic_name = params['ic_name']

slices_dir = run_name+'_slices'
scalars_dir = run_name+'_scalars'
state_dir = run_name+'_state'

dealias = 3/2
timestepper = d3.SBDF2
dtype = np.float64
mesh = None

# Bases
coords = d3.SphericalCoordinates('phi', 'theta', 'r')
dist = d3.Distributor(coords, dtype=dtype, mesh=mesh)
ball = d3.BallBasis(coords, shape=(Nphi, Ntheta, Nr), radius=Radius, dealias=dealias, dtype=dtype)
sphere = ball.surface

# Substitutions
phi, theta, r = dist.local_grids(ball)
lift = lambda A: d3.Lift(A, ball, -1)

# Fields
u = dist.VectorField(coords, name='u',bases=ball)
p = dist.Field(name='p', bases=ball)
tau_p = dist.Field(name='tau_p')
tau_u = dist.VectorField(coords, name='tau u', bases=sphere)
sb_rot = dist.VectorField(coords, name='sb_rot',bases=ball)
sb_rot['g'][0] = (1-ell)*r*np.sin(theta)   # background solid body rotation, np-slip boundary condition, Poincare flow

ex = dist.VectorField(coords, name='ex',bases=ball)
ex['g'][0] = -np.sin(phi)    # phi component of unit vector in x-direction
ex['g'][1] = np.cos(theta)*np.cos(phi)   # theta component of unit vector in x-direction
ex['g'][2] = np.sin(theta)*np.cos(phi)   # r component of unit vector in x-direction

ey = dist.VectorField(coords, name='ey',bases=ball)
ey['g'][0] = np.cos(phi)    # phi component of unit vector in y-direction
ey['g'][1] = np.cos(theta)*np.sin(phi)   # theta component of unit vector in y-direction
ey['g'][2] = np.sin(theta)*np.sin(phi)   # r component of unit vector in y-direction

ez = dist.VectorField(coords, name='ez',bases=ball)
ez['g'][0] = 0    # phi component of unit vector in z-direction
ez['g'][1] = -np.sin(theta)   # theta component of unit vector in x-direction
ez['g'][2] = np.cos(theta)   # r component of unit vector in x-direction


# Problem
problem = d3.IVP([p, u, tau_p, tau_u], namespace=locals())
problem.add_equation("div(u) + tau_p = 0")

# correct equation, but with very small time steps
problem.add_equation("dt(u) - Ek*lap(u) + grad(p) + lift(tau_u) = - u@grad(u) + ell*(2-ell)*ex*(ex@grad(p)) - Ek*ell*(2-ell)*ex@grad(ex@grad(u))")
problem.add_equation("u(r=Radius) = sb_rot(r=Radius)") # Boundary condition
problem.add_equation("integ(p) = 0")  # Pressure gauge

# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

# Initial conditions
if ic_name == 'none':
   u['g'][0]= (1-ell)*r*np.sin(theta) 
   sb_rot['g'][0] = (1-ell)*r*np.sin(theta)
   file_handler_mode = 'overwrite'
else:
    ic_file=ic_name+'_state/'+ic_name+'_state_s1.h5'
    write, initial_time = solver.load_state(ic_file)
    file_handler_mode = 'append'
    print('starting at t = ',solver.sim_time)
    stop_sim_time += solver.sim_time
    solver.stop_sim_time = stop_sim_time

# Analysis
scalars = solver.evaluator.add_file_handler(scalars_dir, iter=50, mode=file_handler_mode)
scalars.add_task(d3.integ( ((u@u)/(1-ell)) + ( ell*(2-ell)*(u@ex)**2/((1-ell)**3) ) )*0.5, layout='g', name='E_kin')

term1 = Ek*( d3.curl(u)@d3.curl(u) )/ (volume*(1-ell))
term2 = Ek*ell*(2-ell)* ( d3.grad(u@ex)@d3.grad(u@ex) ) /(volume*(1-ell)**3)
term3 = -Ek*ell*(2-ell)* ( (ex@d3.grad(u@ex))**2 ) / (volume*(1-ell)**3)
term4 = -Ek*ell*(2-ell)*( (ex@d3.grad(u@ez))**2 ) / (volume*(1-ell))
term5 = -Ek*ell*(2-ell)*( (ex@d3.grad(u@ey))**2 ) / (volume*(1-ell))

scalars.add_task(d3.integ(term1+term2+term3+term4+term5), layout='g', name='E_dis')
scalars.add_task(np.sqrt(u@u), layout='g', name='u_3d')
scalars.add_task(np.sqrt(d3.Curl(u)@d3.Curl(u)), layout='g', name='vorticity_3d')


# Analysis
slices = solver.evaluator.add_file_handler(slices_dir, sim_dt=5, max_writes=200, mode=file_handler_mode)
slices.add_task(u, scales=dealias, name='u_3d') #3D
slices.add_task(d3.Curl(u), name='vorticity_3d') # 3D


# Save final state
final_state = solver.evaluator.add_file_handler(state_dir, sim_dt=stop_sim_time-max_timestep,mode='overwrite')
final_state.add_tasks(solver.state)


# CFL       
CFL = d3.CFL(solver, initial_dt=max_timestep, cadence=10, safety=0.5, threshold=0.05, max_dt=max_timestep)
CFL.add_velocity(u)

# Flow properties
flow = d3.GlobalFlowProperty(solver, cadence=10)
flow.add_property(np.sqrt(u@u), name='v_cfl')

# Main loop
try:
    logger.info('Starting main loop')
    while solver.proceed:
        timestep = CFL.compute_timestep()
        solver.step(timestep)
        #if (solver.iteration-1) % 10 == 0:
        if (solver.iteration) % 100 == 0:
            max_v = flow.max('v_cfl')
            logger.info('Iteration=%i, Time=%e, dt=%e, max(u)=%f' %(solver.iteration, solver.sim_time, timestep, max_v))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()
    







