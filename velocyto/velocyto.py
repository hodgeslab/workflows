#!~/usr/bin/python3

infile = 

n_jobs = 4 # number of threads for kNN calculation

import velocyto as vcy
vlm = vcy.VelocytoLoom(infile)

# Different steps of analysis can be carried on by simply calling the
# methods of this VelocytoLoom object. New variables, normalized version
# of the data matrixes and other parameters will be stored as attributes
# of the “VelocytoLoom” object (method calls will not return any value).
# For example normalization and log transformation can be performed by
# calling the normalize method:

vlm.normalize("S", size=True, log=True)
vlm.S_norm  # contains log normalized

# pdf("1_fractions.pdf", useDingbats=F, height=3, width=3)
vlm.plot_fractions()
# dev.off()

# You can save the results of your analysis in a serialized object
# at any time by running:
# vlm.dump_hdf5("my_velocyto_analysis")

# In another session you can reload the vlm object by running:
# load_velocyto_hdf5("my_velocyto_analysis.hdf5")

# remove the cells with extremely low unspliced detection
vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.5))

# make velocyto aware of the clusters annotation, if we have some
vlm.set_clusters(vlm.ca["ClusterName"])

# Now using the clustering annotation select the genes that are expressed
# above a threshold of total number of molecules in any of the clusters:
vlm.score_detection_levels(min_expr_counts=40, min_cells_express=30)
vlm.filter_genes(by_detection_levels=True)

# perform feature selection:
vlm.score_cv_vs_mean(3000, plot=True, max_expr_avg=35)
vlm.filter_genes(by_cv_vs_mean=True)

# normalize our data by size (total molecule count):
vlm._normalize_S(relative_size=vlm.S.sum(0),
             target_size=vlm.S.sum(0).mean())
vlm._normalize_U(relative_size=vlm.U.sum(0),
             target_size=vlm.U.sum(0).mean())

# For the preparation of the gamma fit we smooth the data using a kNN
# neighbors pooling approach. kNN neighbors can be calculated directly
# in gene expression space or reduced PCA space, using either
# correlation distance or euclidean distance. One example of set of
# parameters is provided below:
vlm.perform_PCA()
vlm.knn_imputation(n_pca_dims=20, k=500, balanced=True, b_sight=3000, b_maxl=1500, n_jobs)

# fit gamma to every gene that survived the filtering step:
vlm.fit_gammas()

# visualize by calling plot_phase_portraits and listing the gene names:
# pdf('2_phase_portraits.pdf', useDingbats=F, height=3, width=3)
vlm.plot_phase_portraits(["AURKA", "NDRG1"])
# dev.off()

# calculate velocity and extrapolate the future state of the cells:
vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1.)

# an alternative extrapolation can be performed using the constant unspliced assumption (for more information consult our preprint)
# vlm.calculate_shift(assumption="constant_unspliced", delta_t=10)
# vlm.extrapolate_cell_at_t(delta_t=1.)

# Projection of velocity onto embeddings. Uuse scikit-learn TSNE
# implementation and make it available as ts attribute as following:
from sklearn.manifold import TSNE
bh_tsne = TSNE()
vlm.ts = bh_tsne.fit_transform(vlm.pcs[:, :25])

# Now we can project on vlm.ts by calling estimate_transition_prob. For
# big datasets this code can take long time to run! We suggest to run it on
# multicore machines (since the implementation is fully multithreaded)
vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1,
                             n_neighbors=3500, knn_random=True, sampled_fraction=0.5)
vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)

# In case of very big dataset visualizations a good way to summarize the
# velocity is to visualize it as velocity field calculated on a grid:

# pdf('3_velocity_field.pdf', useDingbats=F, height=3, width=3)
vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=300)
plt.figure(None,(20,10))
vlm.plot_grid_arrows(quiver_scale=0.6,
                    scatter_kwargs_dict={"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True}, min_mass=24, angles='xy', scale_units='xy',
                    headaxislength=2.75, headlength=5, headwidth=4.8, minlength=1.5,
                    plot_random=True, scale_type="absolute")
# dev.off()
