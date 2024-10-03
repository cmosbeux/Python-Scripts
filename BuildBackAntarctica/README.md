
# Antarctic Geometry Reconstruction

This script reconstructs the ice sheet as it was back in time based on current geometry and dhdt trend. The grounding line position and flotation criterion are based on the hydrostatic equilibrium hypothesis. The reconstruction for a given year is based on the hypothesis that bedmachine geometry is representative of the year 2020.

## Requirements

- **Topography file**: Contains ice thickness and bedrock geometry (e.g., `BedMachineAntarctica-v03.nc`).
- **dhdt files**: For the shelves and the ice sheet, using dhdt observations from IceSat1-2 missions over 2003-2019.

## Potential Requirements

- **Path to data**: Paths are hard-coded for the author's system. You should define the filenames in the different functions.
- **Plotting**: If you want to make some plots, you will need the `Antarctica_Background` module and add it to your paths (available on the author's GitHub).
- **NetCDF handling**: The script uses a custom `format_reading` module for handling NetCDF files (also available on the author's GitHub).
- **Density adjustments**: You might want to adjust water and ice density depending on the densities you will use in your modeling (`ri` and `rw`).

## User Keywords

- `plot`: Boolean to enable/disable plotting.
- `ri`: Ice density (default: 0.917).
- `rw`: Water density (default: 1.028).
- `z_sl`: Sea level (default: 0).
- `targetyear`: The year for reconstruction (default: 1850).
- `resolution`: Resolution in km (default: 1).
- `power_rate`: Rate for long periods, using dhdt rates that decrease with time (default: 1).

## Usage

To run the script, ensure you have the required files and modules, then execute:

```bash
python Antarctic_Geometry.py