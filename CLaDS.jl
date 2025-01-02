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
tree = load_tree("Datasets/Dataset_S12_species_level_phylo.tre")

n_tips(tree)


tip_labels(tree)

# Initialize vector to store sampling fractions
f = a = open("Outputs/CLaDS_sampling_fraction_vector.txt") do f
    readlines(f) |> (s->parse.(Float64, s))
end

# fit ClaDS
output = infer_ClaDS(tree, f =f, print_state = 100, end_it=500000, ltt_steps=1000, thin=1, n_trees=100)

@save "E:/Dropbox/MacroEvoEco_postdoc_Research/projects/Grevilleoideae_biogeography/output/analyses/ClaDS/clads_result_f.jld2" output

save_ClaDS_in_R(output, "E:/Dropbox/MacroEvoEco_postdoc_Research/projects/Grevilleoideae_biogeography/output/analyses/ClaDS/clads_result_f.RData")