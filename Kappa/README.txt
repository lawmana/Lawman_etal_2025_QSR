This code helps reproduce the main figures in "Mechanisms driving the global tropical response to a weakened AMOC during Heinrich Stadial 1" by Lawman et al. publised in Quaternary Science Reviews 

Follow these steps to performing the model-proxy evaluation and then plot figures 1 and 5 of the study:

1. Copy the Lawman_kappa_model.tar.gz file to a working directory

2. Uncompress the .tar file as follows:

	tar -xvf Lawman_kappa_model.tar.gz

3. cd Lawman_kappa_model and run matlab

4. Run kappa_proxy_model_2ensembles.m

This script will perform the proxy-model evaluation reading the synthesis from DiNezio_etal_HS1_database_2025.xlsx and the model data from hose_data_new.mat included in the tar file.
Several flags allow to configure the evaluation for diferent types of proxies and the original interpretation of the records (orig_interp_drywet flag).

5. Generate figure 5

	plot_Fig1b_Lawman_paper
	plot_Fig5_Lawman_paper

