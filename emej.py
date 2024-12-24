# General imports
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple

# Tudat imports
import tudatpy
from tudatpy.trajectory_design import transfer_trajectory
from tudatpy import constants
from tudatpy.numerical_simulation import environment_setup
from tudatpy.util import result2array
from tudatpy.astro.time_conversion import DateTime

# Pygmo imports
import pygmo as pg


# 
def convert_trajectory_parameters(
    transfer_trajectory_object: tudatpy.kernel.trajectory_design.transfer_trajectory.TransferTrajectory,
    trajectory_parameters: List[float],
) -> Tuple[List[float], List[List[float]], List[List[float]]]:
     # Extract node times
    node_times = trajectory_parameters[:4]

    # 提取leg_free_parameters参数
    number_of_leg_free_params = 3  # 有3个leg_free_parameters参数
    leg_free_parameters_start_index = 4
    leg_free_parameters_end_index = leg_free_parameters_start_index + number_of_leg_free_params
    #leg_free_parameters = [trajectory_parameters[leg_free_parameters_start_index:leg_free_parameters_end_index]]
    leg_free_parameters = [[param] for param in trajectory_parameters[leg_free_parameters_start_index:leg_free_parameters_end_index]]

    # 提取node_free_parameters参数
    node_free_parameters = []
    number_of_nodes_minus_one = 4 - 1  # 每个节点（除了最后一个）有3个参数
    node_free_parameters_start_index = leg_free_parameters_end_index
    for i in range(number_of_nodes_minus_one):
        node_free_parameters.append(
            trajectory_parameters[node_free_parameters_start_index + i*3: node_free_parameters_start_index + (i+1)*3]
        )
    # 最后一个节点没有参数，所以添加一个空列表
    node_free_parameters.append([])

    return node_times, leg_free_parameters, node_free_parameters



###########################################################################
# Define transfer trajectory properties
###########################################################################

# Define the central body
central_body = "Sun"

# Define order of bodies (nodes)
transfer_body_order = ["Earth", "Mars", "Earth", "Jupiter"]

# Define departure orbit
departure_semi_major_axis = np.inf
departure_eccentricity = 0

# Define insertion orbit
arrival_semi_major_axis = 1.0895e8 / 0.02
arrival_eccentricity = 0.98

# Create simplified system of bodies
bodies = environment_setup.create_simplified_system_of_bodies()

# Define the trajectory settings for both the legs and at the nodes
transfer_leg_settings, transfer_node_settings = (
    transfer_trajectory.mga_settings_dsm_velocity_based_legs(
        transfer_body_order,
        departure_orbit=(departure_semi_major_axis, departure_eccentricity),
        arrival_orbit=(arrival_semi_major_axis, arrival_eccentricity),
    )
)

# Create the transfer calculation object
transfer_trajectory_object = transfer_trajectory.create_transfer_trajectory(
    bodies,
    transfer_leg_settings,
    transfer_node_settings,
    transfer_body_order,
    central_body,
)

'''
#测试用
trajectory_parameters = [101186930.795904, 135746930.77128, 151158293.62928638,177011172.8831232, 192607302.0999168,
                         0.234594654679, 0.0964769387134, 0.829948744508, 0.317174785637,
                         1408.99421278, 2.3871484244798613, 0.00399193000622855,
                         10931681.13583052, 1.35077257078, 0.0,
                         19397401.558618437, 1.09554368115, 0.0,
                         7015800.056827979, 1.34317576594, 0.0
                         ]


node_times, leg_free_parameters, node_free_parameters = convert_trajectory_parameters(
    transfer_trajectory_object,
    trajectory_parameters
)


print(node_times, leg_free_parameters, node_free_parameters)

'''






#print(node_times, leg_free_parameters, node_free_parameters)

#[101186930.795904, 135746930.77128, 151158293.62928638, 177011172.8831232, 192607302.0999168] [0.234594654679, 0.0964769387134, 0.829948744508, 0.317174785637] [[1408.99421278, 2.3871484244798613, 0.00399193000622855], [10931681.13583052, 1.35077257078, 0.0], [19397401.558618437, 1.09554368115, 0.0], [7015800.056827979, 1.34317576594, 0.0], []]


# Evaluate the transfer with the given parameters

#transfer_trajectory_object.evaluate( node_times, leg_free_parameters, node_free_parameters)


###########################################################################
# CREATE PROBLEM CLASS ####################################################
###########################################################################

class TransferTrajectoryProblem:

    def __init__(
        self,
        transfer_trajectory_object: tudatpy.trajectory_design.transfer_trajectory.TransferTrajectory,
        departure_date_lb: float,  # Lower bound on departure date
        departure_date_ub: float,  # Upper bound on departure date
        legs_tof_lb: np.ndarray,  # Lower bounds of each leg's time of flight
        legs_tof_ub: np.ndarray,  # Upper bounds of each leg's time of flight
        yita_lb: np.ndarray,      # Lower bounds for yita parameters
        yita_ub: np.ndarray,      # Upper bounds for yita parameters
        sudu_lb: np.ndarray,      # Lower bounds for speed parameters
        sudu_ub: np.ndarray,      # Upper bounds for speed parameters
        u_lb: np.ndarray,         # Lower bounds for u parameters
        u_ub: np.ndarray,         # Upper bounds for u parameters
        v_lb: np.ndarray,         # Lower bounds for v parameters
        v_ub: np.ndarray,         # Upper bounds for v parameters
        pe_lb: np.ndarray,        # Lower bounds for pe parameters
        pe_ub: np.ndarray,        # Upper bounds for pe parameters
        jiao_lb: np.ndarray,      # Lower bounds for jiao parameters
        jiao_ub: np.ndarray,       # Upper bounds for jiao parameters
    ):
        self.departure_date_lb = departure_date_lb
        self.departure_date_ub = departure_date_ub
        self.legs_tof_lb = legs_tof_lb
        self.legs_tof_ub = legs_tof_ub
        self.yita_lb = yita_lb   #dsm点火时刻系数
        self.yita_ub = yita_ub
        self.sudu_lb = sudu_lb
        self.sudu_ub = sudu_ub
        self.u_lb = u_lb   #水平角
        self.u_ub = u_ub   #俯仰角
        self.v_lb = v_lb
        self.v_ub = v_ub
        self.pe_lb = pe_lb   #飞跃高度系数
        self.pe_ub = pe_ub
        self.jiao_lb = jiao_lb  #bplane角
        self.jiao_ub = jiao_ub

        self.transfer_trajectory_function = lambda:transfer_trajectory_object

    def get_bounds(self) -> tuple:
        transfer_trajectory_obj = self.transfer_trajectory_function()

        lower_bound = [self.departure_date_lb, self.departure_date_lb + self.legs_tof_lb[0], self.departure_date_lb + self.legs_tof_lb[0] + self.legs_tof_lb[1], self.departure_date_lb + self.legs_tof_lb[0] + self.legs_tof_lb[1] + self.legs_tof_lb[2],
                       self.yita_lb, self.yita_lb, self.yita_lb,
                       self.sudu_lb, self.u_lb, self.v_lb, self.pe_lb * 3.39e6, self.jiao_lb, 0.0, self.pe_lb * 6.378e6, self.jiao_lb, 0.0]
        upper_bound = [self.departure_date_ub, self.departure_date_ub + self.legs_tof_ub[0], self.departure_date_ub + self.legs_tof_ub[0] + self.legs_tof_ub[1], self.departure_date_ub + self.legs_tof_ub[0] + self.legs_tof_ub[1] + self.legs_tof_ub[2],
                       self.yita_ub, self.yita_ub, self.yita_ub,
                       self.sudu_ub, self.u_ub, self.v_ub, self.pe_ub * 3.39e6, self.jiao_ub, 0.0, self.pe_ub * 6.378e6, self.jiao_ub, 0.0]

        #整个问题的难点在这里，这里定义trajectory_parameters里每一个元素的上下界，6.052e6是venus半径，6.378e6是地球半径,3.390e6是火星半径

        bounds = (lower_bound, upper_bound)
        return bounds



    def fitness(self, trajectory_parameters: List[float]) -> list:
        """
        Returns delta V of the transfer trajectory object with the given set of trajectory parameters
        """

        # Retrieve transfer trajectory object
        transfer_trajectory = self.transfer_trajectory_function()

        # Convert list of trajectory parameters to appropriate format
        node_times, leg_free_parameters, node_free_parameters = (
            convert_trajectory_parameters(transfer_trajectory, trajectory_parameters)
        )

        # Evaluate trajectory
        try:
            transfer_trajectory.evaluate(
                node_times, leg_free_parameters, node_free_parameters
            )
            delta_v = transfer_trajectory.delta_v

        # If there was some error in the evaluation of the trajectory, use a very large deltaV as penalty
        except:
            delta_v = 1e10

        return [delta_v]
'''
'''
# 
julian_day = constants.JULIAN_DAY
time1 = (1171.64503236 - 0.5) * julian_day
time2 = time1 + 399.999999715 * julian_day
time3 = time2 + 178.372255301 * julian_day
time4 = time3 + 299.223139512 * julian_day


# 
leg1_param = [0.234594654679]
leg2_param = [0.0964769387134]
leg3_param = [0.829948744508]


#
node1_param = [1408.99421278, 0.37992647165 * 2.0 * np.pi, np.arccos(2.0 * 0.498004040298 - 1.0) - np.pi / 2.0]
node2_param = [1.80629232251 * 3.39e6, 1.35077257078, 0.0]
node3_param = [3.04129845698 * 6.378e6, 1.09554368115, 0.0]
node4_param = []

# 
trajectory_parameters = [
    time1, time2, time3, time4,
    *leg1_param, *leg2_param, *leg3_param,
    *node1_param, *node2_param, *node3_param, *node4_param,
]

'''
'''






# Lower and upper bound on departure date
departure_date_lb = DateTime(2024, 10, 13).epoch()
departure_date_ub = DateTime(2024, 10, 14).epoch()
#print(departure_date_lb)
# List of lower and upper on time of flight for each leg
legs_tof_lb = np.zeros(3)
legs_tof_ub = np.zeros(3)
# mars fly-by
legs_tof_lb[0] = 109 * constants.JULIAN_DAY
legs_tof_ub[0] = 146 * constants.JULIAN_DAY
# earth fly-by
legs_tof_lb[1] = 657 * constants.JULIAN_DAY
legs_tof_ub[1] = 693 * constants.JULIAN_DAY
# jupiter arrive
legs_tof_lb[2] = 1205 * constants.JULIAN_DAY
legs_tof_ub[2] = 1241 * constants.JULIAN_DAY


yita_lb = 0.1
yita_ub = 0.9

sudu_lb = 1000
sudu_ub = 5000

u_lb = 0
u_ub = 3.14

v_lb = 0
v_ub = 3.14

pe_lb = 1
pe_ub = 2

jiao_lb = 0
jiao_ub = 3.14


###########################################################################
# Setup optimization
###########################################################################
# Initialize optimization class
optimizer = TransferTrajectoryProblem(
    transfer_trajectory_object,
    departure_date_lb,
    departure_date_ub,
    legs_tof_lb,
    legs_tof_ub,
    yita_lb,    # Lower bounds for yita parameters
    yita_ub,      # Upper bounds for yita parameters
    sudu_lb,      # Lower bounds for speed parameters
    sudu_ub,      # Upper bounds for speed parameters
    u_lb,         # Lower bounds for u parameters
    u_ub,         # Upper bounds for u parameters
    v_lb,         # Lower bounds for v parameters
    v_ub,         # Upper bounds for v parameters
    pe_lb,        # Lower bounds for pe parameters
    pe_ub,        # Upper bounds for pe parameters
    jiao_lb,      # Lower bounds for jiao parameters
    jiao_ub,
)

# Creation of the pygmo problem object
prob = pg.problem(optimizer)

# To print the problem's information: uncomment the next line
# print(prob)

# Define number of generations per evolution
number_of_generations = 50



# Create pygmo algorithm object
algo = pg.algorithm(pg.sade(gen=number_of_generations))             #自适应差分进化算法

# To print the algorithm's information: uncomment the next line
# print(algo)

# Set population size
population_size = 20

# Create population
pop = pg.population(prob, size=population_size)



###########################################################################
# Run optimization
###########################################################################

# Set number of evolutions
number_of_evolutions = 80

# Initialize empty containers
individuals_list = []
fitness_list = []

for i in range(number_of_evolutions):

    pop = algo.evolve(pop)

    # individuals save
    individuals_list.append(pop.champion_x)
    fitness_list.append(pop.champion_f)

print("The optimization has finished")




###########################################################################
# Results post-processing
###########################################################################

# Extract the best individual
print("\n########### CHAMPION INDIVIDUAL ###########\n")
print("Total Delta V [m/s]: ", pop.champion_f[0])
best_decision_variables = pop.champion_x / constants.JULIAN_DAY
print("Departure time w.r.t J2000 [days]: ", best_decision_variables[0])
print("Earth-Mars time of flight [days]: ", best_decision_variables[1]-best_decision_variables[0])
print("Mars-Earth time of flight [days]: ", best_decision_variables[2]-best_decision_variables[1])
print("Earth-Jupiter time of flight [days]: ", best_decision_variables[3]-best_decision_variables[2])

'''
print('Delta V Timming per leg: ')
print(' - between Earth and Mars:',pop.champion_x[4] * (best_decision_variables[1]-best_decision_variables[0]))

#print('guanjun',pop.champion_x)
'''


'''
#有问题，和best decision vector里的数对不上


print('Delta V per leg: ')
for i in range(len(transfer_body_order)-1):
    print(" - between %s and %s: %.3f m/s" % \
        (transfer_body_order[i], transfer_body_order[i+1], transfer_trajectory_object.delta_v_per_leg[i]))
print()

# Print the DeltaV required at each node
print('Delta V per node : ')
for i in range(len(transfer_body_order)):
    print(" - at %s: %.3f m/s" % \
        (transfer_body_order[i], transfer_trajectory_object.delta_v_per_node[i]))


print('total flytime:',transfer_trajectory_object.time_of_flight/julian_day)
'''



# Plot fitness over generations
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(
    np.arange(0, number_of_evolutions),
    np.float_(fitness_list) / 1000,
    label="Function value: Feval",
)
# Plot champion
champion_n = np.argmin(np.array(fitness_list))



ax.scatter(
    champion_n,
    np.min(fitness_list) / 1000,
    marker="x",
    color="r",
    label="All-time champion",
    zorder=10,
)

# Prettify
ax.set_xlim((0, number_of_evolutions))
ax.set_ylim([10, 50])
ax.grid("major")
ax.set_title("Best individual over generations", fontweight="bold")
ax.set_xlabel("Number of generation")
ax.set_ylabel(r"$\Delta V [km/s]$")
ax.legend(loc="upper right")
plt.tight_layout()

plt.show()

# Reevaluate the transfer trajectory using the champion design variables
node_times, leg_free_parameters, node_free_parameters = convert_trajectory_parameters(
    transfer_trajectory_object, pop.champion_x
)
transfer_trajectory_object.evaluate(
    node_times, leg_free_parameters, node_free_parameters
)

# Extract the state history
state_history = transfer_trajectory_object.states_along_trajectory(500)
fly_by_states = np.array([state_history[node_times[i]] for i in range(len(node_times))])
state_history = result2array(state_history)
au = 1.5e11

# Plot the state history
plt.figure(figsize=(8, 5))
plt.plot(state_history[:, 1] / au, state_history[:, 2] / au)
plt.scatter(
    fly_by_states[0, 0] / au,
    fly_by_states[0, 1] / au,
    color="blue",
    label="Earth departure",
)
plt.scatter(
    fly_by_states[1, 0] / au,
    fly_by_states[1, 1] / au,
    color="green",
    label="Mars fly-by",
)

plt.scatter(
    fly_by_states[2, 0] / au,
    fly_by_states[2, 1] / au,
    color="red",
    label="Earth fly-by",
)
plt.scatter(
    fly_by_states[3, 0] / au,
    fly_by_states[3, 1] / au,
    color="brown",
    label="Jupiter arrive",
)



plt.scatter([0], [0], color="orange", label="Sun")
plt.xlabel("x wrt Sun [AU]")
plt.ylabel("y wrt Sun [AU]")
plt.gca().set_aspect("equal", adjustable='box')
plt.legend(bbox_to_anchor=[1, 1])
plt.show()















