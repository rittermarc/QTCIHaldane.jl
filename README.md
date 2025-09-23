# QTCIhaldane

Code used to generate the plots in the paper  
M. K. Ritter, Y. Núñez Fernández, M. Wallerberger, J. von Delft, H. Shinaoka, and X. Waintal, Quantics Tensor Cross Interpolation for High-Resolution Parsimonious Representations of Multivariate Functions, [Phys. Rev. Lett. 132, 056501](https://doi.org/10.1103/PhysRevLett.132.056501) (2024).

## Usage

### Generating data

Upon being called, the script `scripts/evaluatehaldane.jl` will calculate a Chern number for the Haldane model at some user-specified parameters, and save the corresponding QTT (and some associated information) in a JLD2 data file in the directory specified in `datadirectory`. Note that this directory must exist already, otherwise the script will fail.

Set the parameters `δm`, `R`, and `tolerance` in the file, and execute the script. The main functionality is contained in the function call to `QTCIHaldane.evaluatechern_haldane`, which you can also call from your own scripts in the same manner.

There is an additional boolean parameter `evalooserror`, which indicates whether out-of-sample errors should be evaluated in addition to the in-sample errors used to optimize the QTT. By default, it is set to false, since these additional evaluations significantly increase memory and runtime. You need at least one file with out-of-sample errors to plot these in the second step.

### Plotting data

To plot the data, use the jupyter notebook at `notebook/haldane_plots.ipynb`.

You may have to change the following parameters in cell 2:
- `datadir`: the directory that contains your JLD2 files.
- `nqs`: a list of all `R` values you wish to plot data for.
- `deltams`: a list of all `δm` you wish to plot data for.

Also, in cell 4:
- `datadir_oos`: the directory that contains files with out-of-sample errors.
- Note that a specific file (`R = 20, δm = 1e-5`) will be loaded, change that filename if you need a different one.

Then execute the notebook. It will generate a plot and save it as pdf; by default, the location is `./haldane-6panel.pdf`.

If you used older files, e.g. the ones I originally used for the paper, you may see warnings such as
```
Warning: type TensorCrossInterpolation.TensorCI{Float64} does not exist in workspace; reconstructing
```
These warnings are a consequence of some changes in the TensorCrossInterpolation library, such that the `TensorCI` type is not present any more. The warnings can be safely ignored.
