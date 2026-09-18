.PHONY: figures clean

# The plot scripts write into figures/. Each rebuilds from the aggregated results in
# results*/ and needs none of the primary data; see README.md for the full mapping.
figures:
	Rscript code/sim_rerun_plots.R      # Figures 1, 3 and S1
	Rscript code/fig2_pca.R             # Figure 2
	Rscript code/fig4_prediction.R      # Figure 4
	Rscript code/fig5_real_metrics.R    # Figure 5
	Rscript code/supp2_noqn.R           # Figures S2 and S3
	Rscript code/fig_hist_plot.R        # Figure S4
	Rscript code/supp5_cluster_tile.R   # Figure S5

clean:
	rm -f Rplots.pdf
