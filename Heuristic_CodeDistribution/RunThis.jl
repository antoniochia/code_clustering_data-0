
# set your working directory here

cd("D:\\Trento-2015\\GraphPartitioning\\CodeDistribution")

# Include all subroutine required

include("Set.jl")
include("randomset.jl")
include("CliqueNetHeur.jl")

# read input data

G1=readdlm("G20_1.csv", ';')   # data about connections
C1=readdlm("C20_1.csv", ';')   # data about costs

# set the seed of the random numbers (for reproduceable tests)

srand = 1234

# preliminary operations:

C2 = C1
G2 = makeGsym(int32(rmLastCol(G1)))
N2 = IncToAdj(G2)
numstart = 10*size(G2)[1]

# the two heuristics:

sol1 = CNlocHeur(C2, N2; version = "RR", maxstart = numstart)
sol2 = CNlocHeur(C2, N2; version = "VNS", maxstart = numstart)

# for comments, remarks, further requests and explanation e-mail me at:
# stefano.benati@unitn.it

