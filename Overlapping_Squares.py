import sympy as sp
import matplotlib.pyplot as plt
import numpy as np
from itertools import product

# Define the golden ratio
#phi = (1 + sp.sqrt(5)) / 2
phi = 1.61803398874989484820

# Define the contraction factor
#c = 1 / phi
#c = 3/4
c = 0.70710678118

# Define the variables
x, y = sp.symbols('x y')

# Define the IFS functions
def f_1(x, y):
    return c * x, c * y

def f_2(x, y):
    return c * x, c * y + (1 - c)

def f_3(x, y):
    return c * x + (1 - c), c * y

def f_4(x, y):
    return c * x + (1 - c), c * y + (1 - c)

def generate_addresses(depth):
    addresses = []

    # Generate all possible combinations of digits of length 'depth'
    address_combinations = product(range(1, 5), repeat=depth)

    # Convert combinations to strings and add to the addresses list
    addresses.extend([''.join(map(str, addr)) for addr in address_combinations])

    return addresses

def apply_address_function(j, x, y):
    # Define the IFS functions
    ifs_functions = [f_1, f_2, f_3, f_4]

    # Convert address string to a list of integers
    address = list(map(int, list(j)))

    # Apply the IFS transformations in reverse order
    result = (x, y)
    for idx in reversed(address):
        result = ifs_functions[idx - 1](*result)

    return result

def calculate_tile_corners(j):
    # Define the corners of the unit square
    corners = [(0, 0), (0, 1), (1, 1), (1, 0)]

    # Calculate the corners of the tile for the given address
    tile_corners = [apply_address_function(j, x, y) for x, y in corners]

    return tile_corners

def draw_fractal(depth, filename):
    # Set up the plot
    fig, ax = plt.subplots()
    ax.set_aspect('equal', 'box')
    ax.set_axis_off()

    # Generate addresses in reverse lexicographical order
    addresses = generate_addresses(depth)

    # Draw each tile with random color
    for address in addresses:
        print(address)
        tile_corners = calculate_tile_corners(address)
        # Convert tile corners to numpy array for plotting
        tile_corners_np = np.array(tile_corners)
        print(tile_corners_np)
        # Generate a random color
        color = np.random.rand(3,)

        # Plot the tile
        ax.fill(tile_corners_np[:, 0], tile_corners_np[:, 1], color=color, edgecolor='black', closed=True)
    # Save the figure
    if filename:
        plt.savefig(filename, bbox_inches='tight', pad_inches=0)
    # Show the plot
    plt.show()

# Example: draw the fractal with depth 3
depth = 5
#print(np.array(calculate_tile_corners(generate_addresses(3)[0]))[:,1])
draw_fractal(depth, filename="square-sqrt2-5.png")
