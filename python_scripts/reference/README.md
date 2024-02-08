# Cool & Reference Code


## Plotting
1. - rmsf_gridmap: Grid-based RMSF heatmaps used to show individual residues of beta sheets and lamellar structures. (-> DataFrame)
2. - Heatmap & Annotat_Heatmap: borrowed functions that are used in rmsf_gridmap 
3. - bubble_plot: Seaborn-based scatterplot with (scalable!) sizes proportional to certain weights (-> csv).
4. - old_cv_contour: Uses scipy interpolate on the data.
5. - ddg_scatter: Scatter graph used for molecular weight vs. ddG values. Uses numpy polynomial for basic linear plotting. 

## Analysis
1. - calculate_rmsd & calculate_rgyr: original complete function for calculating and saving the RMSD. Uses the `list(product())` as a list of descriptors, that it then pops and uses as the colums for saving to DataFrame. With usage.   
