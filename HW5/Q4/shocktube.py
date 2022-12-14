import ulula.setups.sedov_taylor as setup_sedov
import ulula.run as ulula_run
import ulula.simulation as ulula_sim
import ulula.plots as ulula_plt
import matplotlib.pyplot as plt
import ulula.setups.shocktube as shock
import numpy as np

# setup the Sod shocktube problem in x-direction
#setup = setup_sedov.SetupSedov()
setup = shock.SetupSodX()

# specify the hydro schemes
hs = ulula_sim.HydroScheme(reconstruction = 'linear', limiter = 'minmod', riemann='hll', time_integration='hancock', cfl = 0.8)

# run the simulation
sim = ulula_run.run(setup, hydro_scheme=hs, tmax=0.2, nx=200)

# plot the 1D images
q_plot = ['DN','VX','PR']
idir=setup.direction
ulula_plt.plot1d(sim, q_plot=q_plot,idir_func=idir)
plt.savefig("shocktube6.png")

# TODO: plot the 2D profiles for density, pressure, and total velocity