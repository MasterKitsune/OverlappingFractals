import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import itertools
from scipy.optimize import fsolve

# Equations for the supergolden ratio and plastic number
def supergolden_equation(x):
    return x**3 - x - 1

def plastic_equation(x):
    return x**3 - x - 1

# Constants
GOLDEN_RATIO = (1 + np.sqrt(5)) / 2
SUPERGOLDEN_RATIO = fsolve(supergolden_equation, 1.0)[0]
PLASTIC_NUMBER = fsolve(plastic_equation, 1.0)[0]
TOLERANCE = 1e-2

# Parameters
ITERATION_DEPTH = 7   # Currently can't go past 10 before decreasing height of tiles
CONTRACTION_VALUE = 1/GOLDEN_RATIO  # Must be between 1/2 and 1

# Functions and Inverses of IFS
def f_1(x):
    return CONTRACTION_VALUE * x

def f_2(x):
    return CONTRACTION_VALUE * x + (1 - CONTRACTION_VALUE)

def f_1_inverse(y):
    return y / CONTRACTION_VALUE

def f_2_inverse(y):
    return (y - (1 - CONTRACTION_VALUE)) / CONTRACTION_VALUE

FUNCTIONS = [f_1, f_2]
INV_FUNCTIONS = [f_1_inverse, f_2_inverse]

# Global Variables
TILE_COLOR_DICT = {}
FIG, AX = plt.subplots()

# Functions for generating Attractor
def lexicographic_list(iteration_depth, num_functions):
    """Generate lexicographically ordered combinations for addressing functions."""
    values = range(1, num_functions + 1)
    return list(itertools.product(values, repeat=iteration_depth))

def apply_ifs_functions(address, point, functions, iteration_depth):
    """Apply iterated function system functions."""
    result = point
    for i in reversed(range(iteration_depth)):
        func_index = address[i] - 1  # Adjust address to 0-based index
        result = functions[func_index](result)
    return result

def apply_ifs_inverse_functions(address, point, inverse_functions, iteration_depth):
    """Apply inverse iterated function system functions."""
    result = point
    for i in range(iteration_depth):
        func_index = address[i] - 1  # Adjust address to 0-based index
        result = inverse_functions[func_index](result)
    return result

def generate_random_color():
    """Generate a random RGB color."""
    return np.random.rand(3,)

def calculate_top_tile(depth, functions, inv_functions, num_functions):
    """Calculate and generate distinct top tiles for each iteration depth."""
    tiles = []

    for i, address in enumerate(lexicographic_list(depth, num_functions)):
        if not tiles:
            L_top = 0  # Set initial value for the first iteration
        else:
            L_top = max([t[1] for t in tiles])
        R_top = apply_ifs_functions(address, 1, functions, depth)

        # Check if the tile is entirely covered by any previous tile
        if any(L <= L_top and R >= R_top for L, R, _ in tiles):
            continue

        # Get the key for the dictionary based on the inverse function
        characterization = apply_ifs_inverse_functions(address, L_top, inv_functions, depth)

        # Check if the tile is equivalent to any previous tile with a tolerance
        equivalent_tile = None
        for key, value in TILE_COLOR_DICT.items():
            if np.allclose(characterization, key, atol=TOLERANCE):  # Adjust tolerance as needed
                equivalent_tile = value
                break

        if equivalent_tile is not None:
            color = equivalent_tile
        else:
            color = generate_random_color()
            TILE_COLOR_DICT[characterization] = color

        tiles.append((L_top, R_top, color))

    return tiles

# Functions for drawing IFS
def draw_top_tiles(top_tiles, height, iteration_depth):
    """Draw rectangles for the calculated top tiles with addresses."""
    for tile in top_tiles:
        color = tile[2]

        # Check if the tile is banned with a tolerance
        if np.isclose(tile[0], tile[1], atol=TOLERANCE) or tile[0] >= tile[1]:
            continue

        rect = patches.Rectangle((tile[0], 1 - height), tile[1] - tile[0], -0.1, linewidth=1,
                                 edgecolor=color, facecolor=color, alpha=0.7)
        AX.add_patch(rect)

        # Add address text inside the rectangle with size decreasing with iteration depth
        if iteration_depth <= 4:
            address_str = ''.join(map(str, lexicographic_list(iteration_depth, len(FUNCTIONS))[top_tiles.index(tile)]))
            text_size = 15 - 2.5 * iteration_depth  # Adjust this value as needed

            plt.text((tile[0] + tile[1]) / 2, 1 - height - 0.05, address_str, fontsize=text_size,
                     ha='center', va='center', color='white' if np.mean(color) < 0.5 else 'black')

    # Add iteration depth label
    plt.text(-0.02, 1 - height - 0.05, f'Iteration {iteration_depth}', fontsize=10, ha='right', va='center')


def draw():
    """Draw the final plot."""
    AX.set_yticks([])
    #AX.set_xlabel('Unit Interval')
    #plt.title('Overlapping Cantor IFS for contraction factor of golden ratio')
    plt.savefig('Cantor-golden.png', bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    for i in range(ITERATION_DEPTH):
        top_tiles = calculate_top_tile(i, FUNCTIONS, INV_FUNCTIONS, len(FUNCTIONS))
        draw_top_tiles(top_tiles, i / ITERATION_DEPTH, i)
    draw()
    print('Number of equivalence Classes:', len(TILE_COLOR_DICT))