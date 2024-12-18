import os
import subprocess
from itertools import product
import re

def run_command(command, capture_output=True):
    """Run a shell command and return its output."""
    result = subprocess.run(command, shell=True, text=True, capture_output=capture_output)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {command}\n{result.stderr}")
    return result.stdout.strip() if capture_output else None

def update_files(tile_x, tile_y, elem_convec_path, main_path):
    """Update tiling parameters for all instances of `tile()` and ensure `nelem=2000`."""
    with open(elem_convec_path, 'r') as file:
        elem_convec_content = file.read()

    # Replace all instances of `tile(...)` with the new parameters
    if "tile(" in elem_convec_content:
        elem_convec_content = re.sub(
            r'tile\(\d+,\d+\)',  # Match current tile(x, y)
            f'tile({tile_x},{tile_y})',  # Replace with the new values
            elem_convec_content
        )
        print(f"Updated tiling parameters to tile({tile_x}, {tile_y}) in {elem_convec_path}.")
    else:
        print("Warning: No tile(...) occurrences found in the file to update.")

    with open(elem_convec_path, 'w') as file:
        file.write(elem_convec_content)

    # Update `nelem` parameter in the main file
    with open(main_path, 'r') as file:
        main_content = file.read()

    main_content = re.sub(
        r'integer\(4\), parameter\s*::\s*nelem\s*=\s*\d+', 
        'integer(4), parameter   :: nelem = 2000', 
        main_content
    )
    print(f"Ensured `nelem` is set to 2000 in {main_path}.")

    with open(main_path, 'w') as file:
        file.write(main_content)

def measure_time():
    """Run the execution script with `time` and parse the real execution time."""
    result = subprocess.run(
        "/usr/bin/time -p ./execution.sh",
        shell=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    print(f"Execution stdout:\n{result.stdout}")
    print(f"Execution stderr:\n{result.stderr}")

    for line in result.stderr.splitlines():
        if line.startswith("real"):
            try:
                # Parse the time and ensure it records three decimal places
                time_value = float(line.split()[1])
                return round(time_value, 3)  # Adjust the number of decimals to 3
            except (ValueError, IndexError):
                raise RuntimeError(f"Failed to parse time from line: {line}")
    raise RuntimeError("Failed to extract execution time from script output.")

# Configuration (hardcoded values)
TILE_X_VALUES = [8 + i for i in range(0, 56, 8)]  # User-defined values
TILE_Y_VALUES = [8 + i for i in range(0, 56, 8)]  # User-defined values
LUMI_BRANCH = "lumi_convec"
TILING_BRANCH = "tiling"
ELEM_CONVEC_FILE = "../src/kernels/sources/elem_convec.f90"
MAIN_FILE = "../src/app/sources/main.f90"
OUTPUT_FILE = "performance_results_cartesian.txt"

def run_experiments(branch):
    """Run experiments for a specific branch."""
    print(f"Switching to branch {branch}...")
    run_command(f"git checkout {branch} --force")

    print(f"Building on branch {branch}...")
    run_command("make clean")
    run_command("make")

    results = []

    for tile_x, tile_y in product(TILE_X_VALUES, TILE_Y_VALUES):
        if tile_x * tile_y > 1024:
            print(f"Skipping tile({tile_x}, {tile_y}) as the product exceeds 1024.")
            continue

        print(f"Running with tile({tile_x}, {tile_y})...")
        update_files(tile_x, tile_y, ELEM_CONVEC_FILE, MAIN_FILE)
        run_command("make clean")
        run_command("make")
        try:
            time = measure_time()
        except Exception as e:
            print(f"Failed to measure time for tile({tile_x}, {tile_y}): {e}")
            time = float('inf')
        results.append((tile_x, tile_y, time))

    return results

def main():
    # Run a single test for lumi_convec branch
    print(f"Switching to branch {LUMI_BRANCH}...")
    run_command(f"git checkout {LUMI_BRANCH} --force")

    print(f"Building on branch {LUMI_BRANCH}...")
    run_command("make clean")
    run_command("make")

    print("Running single test for branch lumi_convec...")
    try:
        time_lumi = measure_time()
    except Exception as e:
        print(f"Failed to measure time for lumi_convec: {e}")
        time_lumi = float('inf')

    # Run experiments on tiling branch
    results_tiling = run_experiments(TILING_BRANCH)

    # Save and process results
    print("Processing results and calculating speedup...")
    with open(OUTPUT_FILE, "w") as file:
        file.write("Tile_X Tile_Y Time_LUMI Time_TILING Speedup\n")
        for tile_x, tile_y, time_tiling in results_tiling:
            speedup = time_lumi / time_tiling if time_tiling > 0 else float('inf')
            file.write(f"{tile_x} {tile_y} {time_lumi:.3f} {time_tiling:.3f} {speedup:.3f}\n")

    print(f"Results saved to {OUTPUT_FILE}")
    with open(OUTPUT_FILE, "r") as file:
        print(file.read())

if __name__ == "__main__":
    main()
