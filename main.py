import matplotlib.pyplot as plt
from math import isclose

# for E  0.005128802728000006 --> 0.005128802729000006 => The wave function has 1 nodes.
# for E  0.005128802729028076, the wave function is well behaved.


x = float(-60)
s = float(0.001)  # s=Δx

k = float(1)  # force_constant
# E = float(input('Energy of Wavefunction: '))

wavefunction = {}
wavefunction[1] = 0
wavefunction[2] = 0.000001

g = {}
dict_for_x = {}


def count_nodes_in_range(wavefunction_values, end_index, start_index=1):
    nodes = 0
    prev_value = wavefunction_values[start_index]

    for next_value in wavefunction_values[start_index + 1: end_index]:
        if prev_value * next_value < 0: nodes = nodes + 1
        prev_value = next_value
    return nodes


a = x
i = 1
while a < -1 * x:
    dict_for_x[i] = a
    a = a + s
    i = i + 1


def calculate_wavefunction(E, dict_for_x):  # dict_for_x is a list now
    for i in dict_for_x:
        g[i] = 2 * 0.00962 - 2 * 0.0000345 * dict_for_x[i] ** 2 + 2 * 0.0000000314 * dict_for_x[i] ** 4 - 2 * E

    for i in range(2, len(dict_for_x)):
        numerator = 2 * wavefunction[i] - wavefunction[i - 1] + (5 * g[i] * wavefunction[i] * s ** 2) / 6 + (
                    g[i - 1] * wavefunction[i - 1] * s ** 2) / 12
        denominator = 1 - (g[i + 1] * s ** 2) / 12
        wavefunction[i + 1] = numerator / denominator

    # end_index = len(wavefunction)
    # num_nodes_in_range = count_nodes_in_range(list(wavefunction.values()), end_index, 1)
    # print(f'The wave function has {num_nodes_in_range} nodes.')
    return wavefunction


def is_well_behaved(wavefunction_values):
    return isclose(wavefunction_values[-1], 0, abs_tol=tolerance)

def find_well_behaved_energy(start_energy, end_energy, energy_step, x_values, user_input_state):
    current_energy = start_energy
    prev_num_nodes = user_input_state
    accuracy_reached = False

    while current_energy <= end_energy + energy_step and not accuracy_reached:
        wavefunction_values = calculate_wavefunction(current_energy, x_values)

        num_nodes_in_range = count_nodes_in_range(list(wavefunction_values.values()), len(wavefunction_values), 1)
        print(f'for E = {current_energy}, The wave function has {num_nodes_in_range} nodes')

        if is_well_behaved(list(wavefunction_values.values())):
            print(f"\nFor E = {current_energy}, THHE WAVEFUNCTION IS WELL BEHAVED WITHIN GIVEN RANGE OF TOLERANCE.")
            end_index = len(wavefunction)
            num_nodes_in_range = count_nodes_in_range(list(wavefunction.values()), end_index, 1)


        if energy_step <= energy_tolerance:
            print(f"Energy accuracy reached: {energy_tolerance}")
            accuracy_reached = True

        if num_nodes_in_range > prev_num_nodes:
            # If the number of nodes changes, print the details and break the loop
            print("\n----------------------------------------------------------------------------------------------------")
            print(f"wave function has energy between {current_energy - energy_step} and {current_energy}.")
            print("----------------------------------------------------------------------------------------------------\n")
            refined_energy = current_energy - energy_step
            new_energy_step = energy_step / 10
            # Refine the energy value with a smaller step for better accuracy
            # Adjust the step size for better accuracy

            if energy_step > energy_tolerance and not accuracy_reached:
                break

        # prev_num_nodes = num_nodes_in_range
        current_energy += energy_step

    if current_energy <= end_energy + energy_step and not accuracy_reached:
        find_well_behaved_energy(refined_energy, refined_energy + energy_step, new_energy_step, x_values, user_input_state)
    return

# User inputs
start_energy = float(input("Enter start energy: "))
end_energy = float(input("Enter end energy: "))
energy_step = float(input("Enter energy step: "))
print("\nThis Method does not work for higher states (More than 10 in case of SHM)")
user_input_state = int(input("Enter user input state: "))
print("\nEnter tolerances and accuracies required")
tolerance = float(input("Enter tolerance for checking the well-behavedness of the wavefunction (suggested: 0.1 to 0.0001): "))
energy_tolerance = float(input("Enter accuracy of the Energy requires (suggested: 0.000001 or more): "))

find_well_behaved_energy(start_energy, end_energy, energy_step, dict_for_x, user_input_state)

# values for double well potential
potential = [0.00962 - 0.0000345 * x ** 2 + 0.0000000314 * x ** 4 for x in dict_for_x.values()]

# Plotting
print("\n\nTrying to plot the wavefunction")
print("Note: Normalisaion takes a lot of time")
plot_speed = int(input("If you want the wavefucntion to be partially scaled down (improper normalisation), enter 0: "))
if plot_speed == 0:
    raw = list(wavefunction.values())
    norm = [500* float(i) / sum(raw) for i in raw]
    plt.plot(list(dict_for_x.values()), norm, linestyle='-', color='blue', label='Ψ(x)', marker='o')
else:
    plt.plot(list(dict_for_x.values()), list(wavefunction.values()), label='Ψ(x)', marker='o')
plt.plot(list(dict_for_x.values()), potential, label='Potential', linestyle='--', color='red')
plt.xlabel('x')
plt.ylabel('Ψ(x)')
plt.title('Ψ(x) vs x')
plt.legend()
plt.grid(which='both', linestyle='-', linewidth='0.5', color='gray')
plt.minorticks_on()
plt.show()