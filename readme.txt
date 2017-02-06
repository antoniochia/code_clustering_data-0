EXACT APPROACHES:



The file FI3+sol_heu.mos is the code of the third incomplete formulation where the solver has been fed with the solution provided by the heuristic approach. 

The file datos20_1.txt is an example of the data file that we have used in our experiments. 

HEURISTIC APPROACHES:


Folder Heuristic_CodeDistribution



JULIA SUBROUTINE CNlocHeur

The Julia code is run on version 0.3.7 (under the JUNO IDE).The subroutine that calculates the optimal partitioning is invoked through:

> sol = CNlocHeur(C, N; version = "xxx", maxstart = yyy)
Compulsory arguments are: 
C = n x n (symmetric) matrix with the cost of pairs (i,j);
N = n x max.degree (adjacency) matrix of a graph G = (V, E).

Optinal arguments are:
version = "RR" or "VNS", depending on which re-start method is implemented;
maxstart = the number of starting solutions, tried to improve.

OUTPUT
sol is a list with: 
1st element: 
0 computation is succesfull	     
>1 there is some mistake somewhere.

2nd element: the value of the objective function

3rd element: the vector v(i), i = 1,...,n with the partition:
 	                v(i) = k stands for "unit i belongs to group k".

4th element: which iteration found the best objective function.	 

EXAMPLE:
The following code is run on Julia 0.3.7.

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

# preliminary operations:C2 = C1

G2 = makeGsym(int32(rmLastCol(G1)))
N2 = IncToAdj(G2)
numstart = 10*size(G2)[1]

# the two heuristics:

sol1 = CNlocHeur(C2, N2; version = "RR", maxstart = numstart)
sol2 = CNlocHeur(C2, N2; version = "VNS", maxstart = numstart)

# for comments, remarks, further requests and explanation e-mail me at:
# stefano.benati@unitn.it 