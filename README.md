# RSC_MEC_curved_surface_navigation

MATLAB code and example data for “Retrosplenial Cortex Exhibits Parallel Directional Representations and Contributes to Medial Entorhinal Spatial Coding in 3D Curved-Surface Navigation”.

## Requirements

MATLAB R2022a was used for code checks. Required toolboxes:

- Statistics and Machine Learning Toolbox.
- Image Processing Toolbox for spatial analyses and Zernike GLM.
- Curve Fitting Toolbox for AIC analysis.

Shared functions, including BNT and CircStat subsets, are included in `functions`.

## General instructions

1. Keep `main`, `functions` and `example_data` in the same repository directory.
2. Open and run an analysis script from `main` in MATLAB. Scripts locate the data and shared functions automatically; analysis settings are defined near the beginning of each script.
3. Run screening before the corresponding downstream analysis:

   - `hd_cylindrical_side_screening.m` before `mlm_cylindrical_side_hd.m`.
   - `hd_platform_screening.m` before `mlm_platform_hd.m`.
   - `position_cylindrical_side_screening.m` before `position_cylindrical_side_1d.m`.

   Dual-axis screening, AIC and Zernike GLM can be run independently.

4. Outputs are saved to a `results` directory created by the scripts. Precomputed results are not included in this repository.

Full-session screening analyses also require `NeuronActivity.mat` in `example_data`.

The `example_data` directory contains:

- `Spilt_behave_calcium_data.mat`: position, head direction, calcium timestamps and calcium events separated into platform and cylindrical-side recordings. `cell_filter_index` contains the original cell indices used for screening.
- `climb_position.mat`: full-session position coordinates and frame indices identifying the platform and cylindrical side.
- `head_direction.mat`: full-session head directions, stored as `cylinder_hd` and `plane_hd`.
- `platform_hd_cell_ids.mat`: a legacy platform HD cell list; current scripts recalculate this list and do not require this file.

`NeuronActivity.mat` contains full-session calcium data; `NeuronActivity.time` contains synchronized timestamps and `Event_filtered_exp2` contains calcium events. Add this file to `example_data` before running full-session analyses.

Times are in seconds, linear positions in centimeters and head directions in degrees. Calcium-event matrices have frames as rows and original cells as columns. In the split data, behavioral and head-direction rows are aligned with calcium frames; the full-session files retain the recording chronology.

Original project code: Copyright (c) 2026 Xuan Zhang and Xin Yuan. Project software is released under the GNU General Public License, version 3; see [LICENSE](LICENSE).

Bundled BNT functions and BNT-derived modifications retain their original notices and GPL terms; see [BNT_LICENSE.txt](BNT_LICENSE.txt). The BNT subset comes from the supplied BNT source tree; upstream project: [cnc-ntnu/BNT](https://bitbucket.org/cnc-ntnu/bnt/).

CircStat functions by Philipp Berens retain their original BSD license in [CircStat_LICENSE.txt](functions/externals/CircStat_LICENSE.txt).

The software license does not cover example data, manuscripts, figures, MATLAB or its toolboxes.
