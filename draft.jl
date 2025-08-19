using CSV
using DataFrames
using Plots
using LaTeXStrings

# Chargement du CSV
df = CSV.read("results_test_3.csv", DataFrame)

plot(df.MTBF_optim, df.cost,
     xlabel="MTBF in optimization",
     ylabel="Simulated cost",
     legend=false,
     marker=:circle,
     label="Cost")

# Ajouter une barre verticale pointillée à true_MTBF
vline!([true_MTBF], linestyle=:dash, color=:red, label="true_MTBF")

annotate!(true_MTBF, -0.1, text(L"\theta^{\mathrm{MTBF}}", :red, 12, :center))


savefig("MTBF_sensitivity.pdf")
display(current())
