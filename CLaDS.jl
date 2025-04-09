# install PANDA
using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("PANDA")
Pkg.add("JLD2")
Pkg.add("RCall")
Pkg.add("RData")

# Library PANDA
using PANDA
using RCall
using RData

# Load the phylo (same as BGB)
tree = load_tree("Dataset_S10.tree")

# Initialize vector to store sampling fractions
# sampling fraction vector generated in BAMM.R script
f = a = open("sampling_fraction_vector.txt") do f
    readlines(f) |> (s->parse.(Float64, s))
end

# fit ClaDS
output = infer_ClaDS(tree, f =f, print_state = 100, end_it=500000, ltt_steps=1000, thin=1, n_trees=100)

# save in R data format
save_ClaDS_in_R(output, "ClaDS_output.RData")