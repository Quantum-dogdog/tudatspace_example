import sys
# Load required standard modules
import os
import numpy as np
from matplotlib import pyplot as plt

# Load required tudatpy modules
from tudatpy import constants
from tudatpy.interface import spice
from tudatpy import numerical_simulation
from tudatpy.numerical_simulation import environment
from tudatpy.numerical_simulation import environment_setup
from tudatpy.numerical_simulation import propagation, propagation_setup
from tudatpy.numerical_simulation import estimation, estimation_setup
from tudatpy.numerical_simulation.estimation_setup import observation
from tudatpy.astro.time_conversion import DateTime
from tudatpy.astro import element_conversion
from tudatpy.util import result2array



# Load spice kernels
path = os.path.abspath('../')
kernels = [path+'/kernels/kernel_noe.bsp']
spice.load_standard_kernels(kernels)

simulation_start_epoch = DateTime(2025,1,1).epoch()
simulation_end_epoch = DateTime(2025,1,10).epoch()

# Define default body settings
bodies_to_create = ["Europa", "Ganymede", "Callisto", "Io", "Jupiter", "Sun"]
global_frame_origin = "Jupiter"
global_frame_orientation = "J2000"
body_settings = environment_setup.get_default_body_settings(bodies_to_create, global_frame_origin,
                                                            global_frame_orientation)

body_settings.add_empty_settings("tianwen4")

# Create system of bodies (in this case only Earth)

bodies = environment_setup.create_system_of_bodies(body_settings)

# Define bodies that are propagated
bodies_to_propagate = ["tianwen4"]

# Define central bodies of propagation
central_bodies = ["Jupiter"]


# Define accelerations acting
acceleration_settings_tianwen4 = dict(
    Europa=[
        propagation_setup.acceleration.spherical_harmonic_gravity(2, 2),
    ],
    Ganymede=[
        propagation_setup.acceleration.spherical_harmonic_gravity(2, 2),
    ],
    Callisto=[
        propagation_setup.acceleration.spherical_harmonic_gravity(2, 2),
    ],
    Io=[
        propagation_setup.acceleration.point_mass_gravity(),
    ],
    Jupiter=[
        propagation_setup.acceleration.spherical_harmonic_gravity(8, 0)
    ],
    Sun=[
        
        propagation_setup.acceleration.point_mass_gravity()
    ]
)

acceleration_settings = {"tianwen4": acceleration_settings_tianwen4}

# Create acceleration models
acceleration_models = propagation_setup.create_acceleration_models(
    bodies, acceleration_settings, bodies_to_propagate, central_bodies
)


initial_state = (280000000,0,0,10000,10000,0)

# Create termination settings
termination_settings = propagation_setup.propagator.time_termination(simulation_end_epoch)


# Create numerical integrator settings
fixed_step_size = 100.0
integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(
    time_step = 100.0,
    coefficient_set = propagation_setup.integrator.rkf_78 )

# Create propagation settings
propagator_settings = propagation_setup.propagator.translational(
    central_bodies,
    acceleration_models,
    bodies_to_propagate,
    initial_state,
    simulation_start_epoch,
    integrator_settings,
    termination_settings
)


# Create simulation object and propagate the dynamics
dynamics_simulator = numerical_simulation.create_dynamics_simulator(
    bodies, propagator_settings
)

# Extract the resulting state history and convert it to an ndarray
states = dynamics_simulator.propagation_results.state_history
states_array = result2array(states)


# # Retrieve the Moon trajectory over vehicle propagation epochs from spice
io_states_from_spice = {
    epoch:spice.get_body_cartesian_state_at_epoch("Io", "Jupiter", "J2000", "None", epoch)
    for epoch in list(states.keys())
}
# # Convert the dictionary to a multi-dimensional array
io_array = result2array(io_states_from_spice)


last_epoch = list(states.keys())[-1]

io_state_at_last_epoch = spice.get_body_cartesian_state_at_epoch("Io", "Jupiter", "J2000", "None", last_epoch)
io_array_at_last_epoch = result2array({last_epoch: io_state_at_last_epoch})

europa_state_at_last_epoch = spice.get_body_cartesian_state_at_epoch("Europa", "Jupiter", "J2000", "None", last_epoch)
europa_array_at_last_epoch = result2array({last_epoch: europa_state_at_last_epoch})

ganymede_state_at_last_epoch = spice.get_body_cartesian_state_at_epoch("Ganymede", "Jupiter", "J2000", "None", last_epoch)
ganymede_array_at_last_epoch = result2array({last_epoch: ganymede_state_at_last_epoch})

callisto_state_at_last_epoch = spice.get_body_cartesian_state_at_epoch("Callisto", "Jupiter", "J2000", "None", last_epoch)
callisto_array_at_last_epoch = result2array({last_epoch: callisto_state_at_last_epoch})

# # Retrieve the Moon trajectory over vehicle propagation epochs from spice
europa_states_from_spice = {
    epoch:spice.get_body_cartesian_state_at_epoch("Europa", "Jupiter", "J2000", "None", epoch)
    for epoch in list(states.keys())
}
# # Convert the dictionary to a multi-dimensional array
europa_array = result2array(europa_states_from_spice)

# # Retrieve the Moon trajectory over vehicle propagation epochs from spice
ganymede_states_from_spice = {
    epoch:spice.get_body_cartesian_state_at_epoch("Ganymede", "Jupiter", "J2000", "None", epoch)
    for epoch in list(states.keys())
}
# # Convert the dictionary to a multi-dimensional array
ganymede_array = result2array(ganymede_states_from_spice)

# # Retrieve the Moon trajectory over vehicle propagation epochs from spice
callisto_states_from_spice = {
    epoch:spice.get_body_cartesian_state_at_epoch("Callisto", "Jupiter", "J2000", "None", epoch)
    for epoch in list(states.keys())
}
# # Convert the dictionary to a multi-dimensional array
callisto_array = result2array(callisto_states_from_spice)

print(
    f"""
Single Earth-Orbiting Satellite Example.
The initial position vector of tw4 is [km]: \n{
    states[simulation_start_epoch][:3] / 1E3}
The initial velocity vector of tw4 is [km/s]: \n{
    states[simulation_start_epoch][3:] / 1E3}
\nAfter {simulation_end_epoch} seconds the position vector of tw4 is [km]: \n{
    states[simulation_end_epoch][:3] / 1E3}
And the velocity vector of tw4 is [km/s]: \n{
    states[simulation_end_epoch][3:] / 1E3}
    """
)


# Define a 3D figure using pyplot
fig = plt.figure(figsize=(10,10), dpi=125)
ax = fig.add_subplot(111, projection='3d')
ax.set_title(f'tw4 trajectory around jupiter')

# Plot the positional state history
ax.plot(states_array[:, 1], states_array[:, 2], states_array[:, 3], label=bodies_to_propagate[0], linestyle='-.')
ax.scatter(0.0, 0.0, 0.0, label="Jupiter", marker='o', color='yellow')
ax.scatter(io_array_at_last_epoch[:,1], io_array_at_last_epoch[:,2], io_array_at_last_epoch[:,3], label="Io", marker='o', color='yellow')
ax.scatter(europa_array_at_last_epoch[:,1], europa_array_at_last_epoch[:,2], europa_array_at_last_epoch[:,3], label="Europa", marker='o', color='yellow')
ax.scatter(ganymede_array_at_last_epoch[:,1], ganymede_array_at_last_epoch[:,2], ganymede_array_at_last_epoch[:,3], label="Ganymede", marker='o', color='yellow')
ax.scatter(callisto_array_at_last_epoch[:,1], callisto_array_at_last_epoch[:,2], callisto_array_at_last_epoch[:,3], label="Callisto", marker='o', color='yellow')
ax.plot(io_array[:,1], io_array[:,2], io_array[:,3], label="Io", linestyle="-", color="grey")
ax.plot(europa_array[:,1], europa_array[:,2], europa_array[:,3], label="Europa", linestyle="-", color="green")
ax.plot(ganymede_array[:,1], ganymede_array[:,2], ganymede_array[:,3], label="Ganymede", linestyle="-", color="blue")
ax.plot(callisto_array[:,1], callisto_array[:,2], callisto_array[:,3], label="Callisto", linestyle="-", color="red")
# Add the legend and labels, then show the plot
ax.legend()
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m]')
plt.show()
