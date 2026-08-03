# MINT Block evaluation

MINT Block is an NP-integration framework for multi-omics data analysis, that deals with the unwanted variation arising from combining data from multiple studies. 

Overview of the MINT Block NP-integration framework. \textbf{a)} MINT Block performs dimensionality reduction analysis for a dataset consisting of observations from $M$ independent studies across $Q$ blocks of variables. \textbf{b)} For each dimension $h \in \{1,...,H\}$ the algorithm computes a (sparse) loading vector $a_h^{(q)}$ and component $t_h^{(q)}$ for each block $q \in \{1,...,Q\}$. After convergence, the input data matrices are deflated and dimension $h + 1$ is computed. \textbf{c)} MINT Block performs multivariate regression or classification analysis, depending on the response variables of interest (continuous or categorial factor), and seeks to maximise the covariance of the components across the blocks.
