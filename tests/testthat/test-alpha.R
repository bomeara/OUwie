context("test-alpha.R")

# Tests issues related to problem demonstrated by Priscilla Lau et al.

test_that("testing likelihood matches", {
	skip_on_cran()


	tree <- structure(
		list(
			edge = structure(
				c(
					6L,
					7L,
					9L,
					9L,
					7L,
					6L,
					8L,
					8L,
					7L,
					9L,
					1L,
					2L,
					3L,
					8L,
					4L,
					5L
				),
				dim = c(8L, 2L)
			),
			edge.length = c(
				0.217684353043458,
				0.679316935884561,
				0.0286901006475091,
				0.0286901006475091,
				0.70800703653207,
				0.675929272806371,
				0.249762116769157,
				0.249762116769157
			),
			tip.label = c("t1", "t2", "t3", "t4", "t5"),
			Nnode = 4L,
			node.label = c(1, 1, 2, 1)
		),
		class = "phylo",
		order = "cladewise"
	)

	trait <- structure(
		list(
			Genus_species = c("t1", "t2", "t3", "t4", "t5"),
			Reg = c(1, 1, 1, 2, 2),
			X = c(0.27768857, 0.08634373, 0.02390178, 1.71202848, 1.91449654)
		),
		class = "data.frame",
		row.names = c(NA, -5L)
	)

	tree$edge.length <- tree$edge.length * 20
	disc_trait <- trait$Reg
	cont_trait <- trait$X
	names(disc_trait) <- names(cont_trait) <- trait$Genus_species
	set.seed(1137)
		
	tree <- invisible(phytools::make.simmap(tree, disc_trait)) 

	# OUMVA
	alpha = c(0.5632459, 0.1726052)
	sigma.sq = c(0.1064417, 0.3461386)
	theta = c(1.678196, 0.4185894)

	names(alpha) <- names(sigma.sq) <- names(theta) <- c("1", "2")

	alpha_full <- alpha
	sigma_full <- sigma.sq
	theta_full <- theta

	OUMVA_ouwie <- OUwie.fixed(
		tree,
		trait,
		model = c("OUMVA"),
		simmap.tree = TRUE,
		scaleHeight = FALSE,
		clade = NULL,
		alpha = alpha,
		sigma.sq = sigma.sq,
		theta = theta,
		algorithm = "invert"
	)
	
	lau_loglikelhood <- sd_logL_vcv(
		tree,
		cont_trait,
		alpha,
		sigma.sq,
		theta
	)
	
	

	expect_equal(lau_loglikelhood, OUMVA_ouwie$loglik[1,1], tolerance=1e-3)
})


test_that("testing vcv oumva", {
	skip_on_cran()

	tree <- structure(
		list(
			edge = structure(
				c(
					6L,
					7L,
					9L,
					9L,
					7L,
					6L,
					8L,
					8L,
					7L,
					9L,
					1L,
					2L,
					3L,
					8L,
					4L,
					5L
				),
				dim = c(8L, 2L)
			),
			edge.length = c(
				0.217684353043458,
				0.679316935884561,
				0.0286901006475091,
				0.0286901006475091,
				0.70800703653207,
				0.675929272806371,
				0.249762116769157,
				0.249762116769157
			),
			tip.label = c("t1", "t2", "t3", "t4", "t5"),
			Nnode = 4L,
			node.label = c(1, 1, 2, 1)
		),
		class = "phylo",
		order = "cladewise"
	)

	trait <- structure(
		list(
			Genus_species = c("t1", "t2", "t3", "t4", "t5"),
			Reg = c(1, 1, 1, 2, 2),
			X = c(0.27768857, 0.08634373, 0.02390178, 1.71202848, 1.91449654)
		),
		class = "data.frame",
		row.names = c(NA, -5L)
	)

	tree$edge.length <- tree$edge.length * 20
	disc_trait <- trait$Reg
	cont_trait <- trait$X
	names(disc_trait) <- names(cont_trait) <- trait$Genus_species
	set.seed(1137)

	tree <- invisible(phytools::make.simmap(tree, disc_trait))

	# OUMVA
	alpha = c(0.5632459, 0.1726052)
	sigma.sq = c(0.1064417, 0.3461386)
	theta = c(1.678196, 0.4185894)

	names(alpha) <- names(sigma.sq) <- names(theta) <- c("1", "2")

	alpha_full <- alpha
	sigma_full <- sigma.sq
	theta_full <- theta

	OUMVA_ouwie <- OUwie.fixed(
		tree,
		trait,
		model = c("OUMVA"),
		simmap.tree = TRUE,
		scaleHeight = FALSE,
		clade = NULL,
		alpha = alpha,
		sigma.sq = sigma.sq,
		theta = theta,
		algorithm = "invert"
	)

	lau_vcv <- vcv.matrix(tree, alpha, sigma.sq)

	root.edge.index <- which(tree$edge[, 1] == ape::Ntip(tree) + 1)

	root.state <- which(
		colnames(tree$mapped.edge) == names(tree$maps[[root.edge.index[2]]][1])
	)

	ouwie_vcv <- OUwie:::varcov.ou.enhanced.tree(
		phy = tree,
		enhanced_tree = OUwie:::create_enhanced_tree_structure(tree),
		Rate.mat = OUMVA_ouwie$solution,
		root.state = root.state,
		root.age = NULL,
		scaleHeight = FALSE
	)

	expect_equal(
		lau_vcv,
		ouwie_vcv,
		tolerance = 1e-3
	)
})



test_that("testing vcv bm1", {
	skip_on_cran()

	tree <- structure(
		list(
			edge = structure(
				c(
					6L,
					7L,
					9L,
					9L,
					7L,
					6L,
					8L,
					8L,
					7L,
					9L,
					1L,
					2L,
					3L,
					8L,
					4L,
					5L
				),
				dim = c(8L, 2L)
			),
			edge.length = c(
				0.217684353043458,
				0.679316935884561,
				0.0286901006475091,
				0.0286901006475091,
				0.70800703653207,
				0.675929272806371,
				0.249762116769157,
				0.249762116769157
			),
			tip.label = c("t1", "t2", "t3", "t4", "t5"),
			Nnode = 4L,
			node.label = c(1, 1, 2, 1)
		),
		class = "phylo",
		order = "cladewise"
	)

	trait <- structure(
		list(
			Genus_species = c("t1", "t2", "t3", "t4", "t5"),
			Reg = c(1, 1, 1, 2, 2),
			X = c(0.27768857, 0.08634373, 0.02390178, 1.71202848, 1.91449654)
		),
		class = "data.frame",
		row.names = c(NA, -5L)
	)

	tree$edge.length <- tree$edge.length * 20
	disc_trait <- trait$Reg
	cont_trait <- trait$X
	names(disc_trait) <- names(cont_trait) <- trait$Genus_species
	set.seed(1137)
	tree <- invisible(phytools::make.simmap(tree, disc_trait))

	# BM1
	alpha = rep(1e-8, 2)
	sigma.sq = c(1, 1)
	theta = c(42, 42)

	names(alpha) <- names(sigma.sq) <- names(theta) <- c("1", "2")

	alpha_full <- alpha
	sigma_full <- sigma.sq
	theta_full <- theta

	lau_vcv <- vcv.matrix(tree, alpha, sigma.sq)

	root.edge.index <- which(tree$edge[, 1] == ape::Ntip(tree) + 1)

	root.state <- which(
		colnames(tree$mapped.edge) == names(tree$maps[[root.edge.index[2]]][1])
	)

	enhanced_tree <- OUwie:::create_enhanced_tree_structure(tree)
	
	Rate.mat <- matrix(0, nrow=2, ncol=2)
	Rate.mat[1,] <- alpha
	Rate.mat[2,] <- sigma.sq
	ouwie_vcv <- OUwie:::varcov.ou.enhanced.tree(
		phy = tree,
		enhanced_tree = enhanced_tree,
		Rate.mat = Rate.mat,
		root.state = root.state,
		root.age = NULL,
		scaleHeight = FALSE
	)

	expect_equal(
		lau_vcv,
		ouwie_vcv,
		tolerance = 1e-3
	)
})


