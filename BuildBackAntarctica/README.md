
# Antarctic Geometry Reconstruction

This script reconstructs the ice sheet as it was back in time based on current geometry and dhdt trend. The grounding line position and flotation criterion are based on the hydrostatic equilibrium hypothesis. The reconstruction for a given year is based on the hypothesis that bedmachine geometry is representative of the year 2020.

## Requirements

- **Topography file**: Contains ice thickness and bedrock geometry (e.g., `BedMachineAntarctica-v03.nc`).
- **dhdt files**: For the shelves and the ice sheet, using dhdt observations from IceSat1-2 missions over 2003-2019.

> **⚠️ Warning:**
> Make sure to have the correct path for the different data files. The paths are hard-coded into the Python script.

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
```

## Functions

Build_IceSat2_ContinuousMap(x_final=None, y_final=None, masked=False)
Reads dHdT from IceSat1-2 over the shelves and builds a continuous map.

### Parameters

- `x_final`: Final x-coordinates (default: None).
- `y_final`: Final y-coordinates (default: None).
- `masked`: Boolean to apply masking (default: False).

### Example

```python
Build_IceSat2_ContinuousMap(x_final=some_x, y_final=some_y, masked=True)
```

## Dependencies

- `numpy`
- `matplotlib`
- `scipy`
- `h5py`
- `Antarctica_Background`(custom module available in MyModules)
- `format_reading`(custom module available in MyModules)
- `FunctionsDataPlot` (custom module in the same folder as `Antarctic_Geometry.py`)

## Author

- **cmosbeux**

## License

This project is licensed under the MIT License - see the LICENSE file for details.