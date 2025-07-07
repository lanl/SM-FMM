rm -r figures
mkdir figures

cp spherical_inclusion/FMM/micro_coarse/sph_initial_micro.png ./figures/
cp spherical_inclusion/FMM/micro_coarse/sph_line_plot_strain_33.png ./figures/
cp spherical_inclusion/FMM/micro_coarse/sph_str_dif.png ./figures/
cp spherical_inclusion/FMM/micro_element/sph_str_dif_nel.png ./figures/

cp px/FMM/vary_box/px_FMMbx_run_time_ord0.png ./figures/
cp px/FMM/vary_box/px_FMMbx_run_time_ord1.png ./figures/
cp px/FMM/vary_box/px_FMMbx_run_time_ord2.png ./figures/
cp px/FMM/vary_box/px_FMMbx_run_time_ord3.png ./figures/
cp px/FMM/vary_box/px_FMMbx_run_time_ord4.png ./figures/

cp px/px_epvm_0020.png ./figures/
cp px/px_run_time.png ./figures/
cp px/px_stress_hyd_0020.png ./figures/
cp px/px_stress_strain.png ./figures/
cp px/px_stressZZ_0001.png ./figures/
cp px/px_svm_0020.png ./figures/

cp px_cracked_hcp/pxc_initial_micro.png ./figures/
cp px_cracked_hcp/pxc_iso_stiff_vs_thick.png ./figures/
cp px_cracked_hcp/pxc_compare_svm.png ./figures/