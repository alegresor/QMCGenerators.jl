var documenterSearchIndex = {"docs":
[{"location":"#Quasi-Monte-Carlo-Generators","page":"Home","title":"Quasi-Monte Carlo Generators","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Julia version of Dirk Nuyens' Magic Point Shop. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the Tutorial for instructions on how to use this package.  ","category":"page"},{"location":"#Quasi-Random-Generators","page":"Home","title":"Quasi-Random Generators","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"QMCGenerators.jl implements quasi-random (low discrepancy) sequence generators. Quasi-random points are carefully coordinated in a dependent manner to more evenly fill 01^s compared to pseudo-random (independent) points. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package implements two flavors of quasi-random sequences: Lattice rules and digital nets, both in base 2. Independent randomizations may be applied to base sequences via random shifts for Lattices and random digital shifts for digital nets.  ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: image)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Notice the gaps and clusters of the above pseudo-random points compared to the more even coverage of quasi-random points. The implemented quasi-random sequences are extensible, meaning you may increase the sample size through powers of two to better fill the space. Above, each plot starts with 2^6 blue points. Doubling the sample size to 2^7 adds the green points. Doubling again to 2^8 samples adds the purple points. ","category":"page"},{"location":"#Quasi-Monte-Carlo","page":"Home","title":"Quasi-Monte Carlo","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Monte Carlo methods approximate the true mean","category":"page"},{"location":"","page":"Home","title":"Home","text":"mu = mathbbEf(X) = int_01^s f(x) mathrmd x","category":"page"},{"location":"","page":"Home","title":"Home","text":"for f 01^s to mathbbR and X sim mathcalU01^s by the sample mean","category":"page"},{"location":"","page":"Home","title":"Home","text":"hatmu = frac1n sum_i=0^n-1 f(X_i)","category":"page"},{"location":"","page":"Home","title":"Home","text":"where X_0X_1dots sim mathcalU01^s. If X_i are chosen to be pseudo-random, then the sample average is a Simple Monte Carlo approximation with error mathcalO(n^-12). If instead we choose X_i to be quasi-random then the sample average is a Quasi-Monte Carlo approximation with error approaching mathcalO(n^-1). ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: image)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The above example approximates mu = mathbbEg(T) where g mathbbR^s to mathbbR is the Keister function","category":"page"},{"location":"","page":"Home","title":"Home","text":"g(t) = pi^s2 cos( lVert t rVert_2)","category":"page"},{"location":"","page":"Home","title":"Home","text":"and T sim mathcalN(0I_s2)  in dimension s=7. To make this Quasi-Monte Carlo compatible, we set ","category":"page"},{"location":"","page":"Home","title":"Home","text":"f(x) = g(Phi^-1(x)sqrt2)","category":"page"},{"location":"","page":"Home","title":"Home","text":"where Phi^-1 is the inverse CDF of mathcalN(01) applied elementwise. This transform satisfies","category":"page"},{"location":"","page":"Home","title":"Home","text":"mu = mathbbEf(X) = mathbbEg(T)","category":"page"},{"location":"","page":"Home","title":"Home","text":"where X sim mathcalU01^s as desired. ","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Kuo, F. Y., & Nuyens, D. (2016). Application of quasi-Monte Carlo methods to elliptic PDEs with random diffusion coefficients: a survey of analysis and implementation. Foundations of Computational Mathematics, 16, 1631-1696.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"DocTestSetup = quote\n    using QMCGenerators\n    using CairoMakie\n    using Random: MersenneTwister\n    using Distributions\n    using LinearAlgebra\nend","category":"page"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Pages = [\"tutorial.md\"]\nDepth = 4","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To begin, install this package with ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Pkg.add(\"QMCGenerators\")","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"and then import via","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using QMCGenerators","category":"page"},{"location":"tutorial/#Common-Usage","page":"Tutorial","title":"Common Usage","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Often we just want to generate a single randomized quasi-random sequence. For example, to generate the first 4 points of a 3 dimensional digital net with 1 random digital shift with seed 7","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Next(RandomDigitalShift(DigitalSeqB2G(3),1,7),4)\n# output\n4×3 Matrix{Float64}:\n 0.243795  0.719182  0.810814\n 0.743795  0.219182  0.310814\n 0.993795  0.969182  0.560814\n 0.493795  0.469182  0.0608139","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"A similar API is available for randomly shifted Lattices:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Next(RandomShift(LatticeSeqB2(3),1,7),4)\n# output\n4×3 Matrix{Float64}:\n 0.243795  0.719182  0.810814\n 0.743795  0.219182  0.310814\n 0.493795  0.969182  0.0608139\n 0.993795  0.469182  0.560814","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"warning: Warning\nWhile not strictly enforced, sample sizes should be powers of two to achieve full coverage of 01^s","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"In the following sections we always supply a seed for reproducibility. Supplying a seed requires you also supply the number of randomizations as done above. However, if you do not wish to seed, you can simply supply the number of randomizations. ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"rls = RandomDigitalShift(DigitalSeqB2G(12),2)\nxs = NextR(rls,2^7)\nsize(xs)\n# output\n(2,)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"size(xs[1]) == size(xs[2]) == (2^7,12)\n# output\ntrue","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Moreover, you may use only one randomization without a seed with the simplified API","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"rls = RandomDigitalShift(DigitalSeqB2G(52))\nx = Next(rls,2^14)\nsize(x)\n# output\n(16384, 52)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The same API simplifications holds for Lattices. ","category":"page"},{"location":"tutorial/#Structure-and-Functions","page":"Tutorial","title":"Structure and Functions","text":"","category":"section"},{"location":"tutorial/#Unrandomized-Sequences","page":"Tutorial","title":"Unrandomized Sequences","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"warning: Warning\nIt is highly recommended you randomize sequences. The first point of an unrandomized sequence is 0 in 01^s which will be transformed to an infinite values in many functions e.g. those composed with taking the inverse CDF of a mathcalN(01) as done in the section on Quasi-Monte Carlo. ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Let's start by defining a 5 dimensional (unrandomized) digital sequence and generating the first 4 points.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"ds = DigitalSeqB2G(5)\nNext(ds,4)\n# output\n4×5 Matrix{Float64}:\n 0.0   0.0   0.0   0.0   0.0\n 0.5   0.5   0.5   0.5   0.5\n 0.75  0.25  0.25  0.25  0.75\n 0.25  0.75  0.75  0.75  0.25","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"We can then generate the next 4 points in the sequence with ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Next(ds,4)\n# output\n4×5 Matrix{Float64}:\n 0.375  0.375  0.625  0.875  0.375\n 0.875  0.875  0.125  0.375  0.875\n 0.625  0.125  0.875  0.625  0.625\n 0.125  0.625  0.375  0.125  0.125","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To reset the generator use","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Reset!(ds)\nNext(ds,4)\n# output\n4×5 Matrix{Float64}:\n 0.0   0.0   0.0   0.0   0.0\n 0.5   0.5   0.5   0.5   0.5\n 0.75  0.25  0.25  0.25  0.75\n 0.25  0.75  0.75  0.75  0.25","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Let's reset once more before continuing ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Reset!(ds)\n# output","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"These functions can also be applied to Lattices ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"ls = LatticeSeqB2(5)\nNext(ls,4)\n# output\n4×5 Matrix{Float64}:\n 0.0   0.0   0.0   0.0   0.0\n 0.5   0.5   0.5   0.5   0.5\n 0.25  0.25  0.25  0.25  0.25\n 0.75  0.75  0.75  0.75  0.75","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Next(ls,4)\n# output\n4×5 Matrix{Float64}:\n 0.125  0.625  0.125  0.625  0.625\n 0.625  0.125  0.625  0.125  0.125\n 0.375  0.875  0.375  0.875  0.875\n 0.875  0.375  0.875  0.375  0.375","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Reset!(ls)\n# output","category":"page"},{"location":"tutorial/#Independent-Sequence-Randomizations","page":"Tutorial","title":"Independent Sequence Randomizations","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Getting a single randomization was shown in the Common Usage section. Support is also available for multiple independent randomizations. For instance, we can generating 2 independent randomizations with seed 11 and get the next 4 points in each sequence via  ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"rds = RandomDigitalShift(ds,2,11)\nxs = NextR(rds,4)\ntypeof(xs)\n# output\nVector{Matrix{Float64}} (alias for Array{Array{Float64, 2}, 1})","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"xs[1]\n# output\n4×5 Matrix{Float64}:\n 0.58051    0.862848  0.861659  0.919821  0.35264\n 0.0805097  0.362848  0.361659  0.419821  0.85264\n 0.33051    0.612848  0.611659  0.669821  0.60264\n 0.83051    0.112848  0.111659  0.169821  0.10264","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"xs[2]\n# output\n4×5 Matrix{Float64}:\n 0.714139  0.702579  0.0130147  0.901326  0.280776\n 0.214139  0.202579  0.513015   0.401326  0.780776\n 0.464139  0.952579  0.263015   0.651326  0.530776\n 0.964139  0.452579  0.763015   0.151326  0.0307758","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"As with unrandomized sequences, we can get the next 4 points and then reset the generator with ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"xs = NextR(rds,4)\nxs[1]\n# output\n4×5 Matrix{Float64}:\n 0.95551  0.737848  0.486659  0.0448206  0.22764\n 0.45551  0.237848  0.986659  0.544821   0.72764\n 0.20551  0.987848  0.236659  0.294821   0.97764\n 0.70551  0.487848  0.736659  0.794821   0.47764","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"xs[2]\n# output\n4×5 Matrix{Float64}:\n 0.839139   0.827579   0.638015  0.0263257  0.155776\n 0.339139   0.327579   0.138015  0.526326   0.655776\n 0.0891388  0.577579   0.888015  0.276326   0.905776\n 0.589139   0.0775787  0.388015  0.776326   0.405776","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Reset!(rds)\n# output","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Similarly for Lattices ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"rls = RandomShift(ls,2,11)\nxr = NextR(rls,4)\nxr[1]\n# output\n4×5 Matrix{Float64}:\n 0.498434  0.26454    0.676602  0.46979  0.677608\n 0.998434  0.76454    0.176602  0.96979  0.177608\n 0.748434  0.51454    0.926602  0.71979  0.927608\n 0.248434  0.0145404  0.426602  0.21979  0.427608","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"xr[2]\n# output\n4×5 Matrix{Float64}:\n 0.389721  0.719424  0.184079  0.568002   0.105514\n 0.889721  0.219424  0.684079  0.0680019  0.605514\n 0.639721  0.969424  0.434079  0.818002   0.355514\n 0.139721  0.469424  0.934079  0.318002   0.855514","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Reset!(rls)\n# output","category":"page"},{"location":"tutorial/#Advanced-Features","page":"Tutorial","title":"Advanced Features","text":"","category":"section"},{"location":"tutorial/#Alternative-Generating-Matrices-and-Vectors","page":"Tutorial","title":"Alternative Generating Matrices and Vectors","text":"","category":"section"},{"location":"tutorial/#Pregenerated","page":"Tutorial","title":"Pregenerated","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"We include many alternative generating matrices for digital sequences in this directory. To use an alternative, simply supply the relative path ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"ds = DigitalSeqB2G(3,\"sobolmats/sobol_alpha2_Bs64.col\")\nNext(ds,4)\n# output\n4×3 Matrix{Float64}:\n 0.0     0.0     0.0\n 0.75    0.75    0.75\n 0.6875  0.1875  0.9375\n 0.4375  0.9375  0.1875","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Alternative Lattice generating vectors are available in this directory. For Lattices, after supplying the path you also need to pass the m value in the file name","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"ls = LatticeSeqB2(3,\"exod8_base2_m13.txt\",13)\nNext(ls,4)\n# output\n4×3 Matrix{Float64}:\n 0.0   0.0   0.0\n 0.5   0.5   0.5\n 0.25  0.25  0.75\n 0.75  0.75  0.25","category":"page"},{"location":"tutorial/#User-Defined","page":"Tutorial","title":"User Defined","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"One may supply their own generating matrix to construct a base 2 digital sequence, for example","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"m = 5\nC1 = [BigInt(2^i) for i=0:(m-1)]\nC2 = [BigInt(1) for i=1:m]\nfor i in 2:m C2[i] = (C2[i-1] << 1) ⊻ C2[i-1] end\ngenerating_matrix = vcat(C1',C2')\n# output\n2×5 Matrix{BigInt}:\n 1  2  4   8  16\n 1  3  5  15  17","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"ds = DigitalSeqB2G(2,generating_matrix)\nNext(ds,4)\n# output\n4×2 Matrix{Float64}:\n 0.0   0.0\n 0.5   0.5\n 0.75  0.25\n 0.25  0.75","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"For base 2 Lattices, you may supply the generating vector followed by m where 2^m is the maximum number of supported points ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"generating_vector = BigInt[1,433461,315689]\nm = 20\nls = LatticeSeqB2(3,generating_vector,m)\nNext(ls,4)\n# output\n4×3 Matrix{Float64}:\n 0.0   0.0   0.0\n 0.5   0.5   0.5\n 0.25  0.25  0.25\n 0.75  0.75  0.75","category":"page"},{"location":"tutorial/#Linear-Ordering","page":"Tutorial","title":"Linear Ordering","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"By default, digital sequences are generated in Gray code order. One may generate the first 2^m points in linear order via ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"m = 3\nn = 2^m\nds = DigitalSeqB2G(4)\nFirstLinear(ds,m)\n# output\n8×4 Matrix{Float64}:\n 0.0    0.0    0.0    0.0\n 0.5    0.5    0.5    0.5\n 0.25   0.75   0.75   0.75\n 0.75   0.25   0.25   0.25\n 0.125  0.625  0.375  0.125\n 0.625  0.125  0.875  0.625\n 0.375  0.375  0.625  0.875\n 0.875  0.875  0.125  0.375","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Compare to the original ordering ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Next(ds,n)\n# output\n8×4 Matrix{Float64}:\n 0.0    0.0    0.0    0.0\n 0.5    0.5    0.5    0.5\n 0.75   0.25   0.25   0.25\n 0.25   0.75   0.75   0.75\n 0.375  0.375  0.625  0.875\n 0.875  0.875  0.125  0.375\n 0.625  0.125  0.875  0.625\n 0.125  0.625  0.375  0.125","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Similarly, Lattices are by default generated in extensible ordering. A linear ordering is also available ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"m = 3\nn = 2^m \nls = LatticeSeqB2(4)\nFirstLinear(ls,m)\n# output\n8×4 Matrix{Float64}:\n 0.0    0.0    0.0    0.0\n 0.125  0.625  0.125  0.625\n 0.25   0.25   0.25   0.25\n 0.375  0.875  0.375  0.875\n 0.5    0.5    0.5    0.5\n 0.625  0.125  0.625  0.125\n 0.75   0.75   0.75   0.75\n 0.875  0.375  0.875  0.375","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Compare to the original ordering ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Next(ls,n)\n# output\n8×4 Matrix{Float64}:\n 0.0    0.0    0.0    0.0\n 0.5    0.5    0.5    0.5\n 0.25   0.25   0.25   0.25\n 0.75   0.75   0.75   0.75\n 0.125  0.625  0.125  0.625\n 0.625  0.125  0.625  0.125\n 0.375  0.875  0.375  0.875\n 0.875  0.375  0.875  0.375","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Linear order for randomized sequences has expected syntax ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"ds = DigitalSeqB2G(4)\nrds = RandomDigitalShift(ds,1,17)\nFirstLinear(rds,2)\n# output\n4×4 Matrix{Float64}:\n 0.844967   0.66901  0.686087  0.362606\n 0.344967   0.16901  0.186087  0.862606\n 0.594967   0.41901  0.436087  0.612606\n 0.0949675  0.91901  0.936087  0.112606","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"rds = RandomDigitalShift(ds,2,17)\nxs = FirstRLinear(rds,2)\nxs[1]\n# output\n4×4 Matrix{Float64}:\n 0.844967   0.686087  0.936228  0.96835\n 0.344967   0.186087  0.436228  0.46835\n 0.594967   0.436087  0.186228  0.21835\n 0.0949675  0.936087  0.686228  0.71835","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"xs[2]\n# output\n4×4 Matrix{Float64}:\n 0.66901  0.362606  0.532599   0.347617\n 0.16901  0.862606  0.0325994  0.847617\n 0.91901  0.612606  0.282599   0.597617\n 0.41901  0.112606  0.782599   0.0976169","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The same functions are available for randomly shifted lattices. ","category":"page"},{"location":"tutorial/#Binary-Functions-for-Digital-Sequences","page":"Tutorial","title":"Binary Functions for Digital Sequences","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"For digital sequences, we sometimes want the binary representation of points. We can get the binary representations as integers and then convert them to their floating point values as follows ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"ds = DigitalSeqB2G(4)\nxb = NextBinary(ds,4)\n# output\n4×4 Matrix{BigInt}:\n          0           0           0           0\n 2147483648  2147483648  2147483648  2147483648\n 3221225472  1073741824  1073741824  1073741824\n 1073741824  3221225472  3221225472  3221225472","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"BinaryToFloat64(xb,ds)\n# output\n4×4 Matrix{Float64}:\n 0.0   0.0   0.0   0.0\n 0.5   0.5   0.5   0.5\n 0.75  0.25  0.25  0.25\n 0.25  0.75  0.75  0.75","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"This is also compatible with randomized digital sequences ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Reset!(ds)\nrds_single = RandomDigitalShift(ds,1,11)\nxb = NextBinary(rds_single,4)\n# output\n4×4 Matrix{BigInt}:\n 5228766587057444  6432390138047189  7771840027009430  6328265944608255\n  725166959686948  1928790510676693  3268240399638934  1824666317237759\n 2976966773372196  8684189951732437  5520040213324182  8580065758293503\n 7480566400742692  4180590324361941  1016440585953686  4076466130923007","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"BinaryToFloat64(xb,rds_single)\n# output\n4×4 Matrix{Float64}:\n 0.58051    0.714139  0.862848  0.702579\n 0.0805097  0.214139  0.362848  0.202579\n 0.33051    0.964139  0.612848  0.952579\n 0.83051    0.464139  0.112848  0.452579","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Reset!(ds)\nrds_multiple = RandomDigitalShift(ds,2,11)\nxbs = NextRBinary(rds_multiple,4)\nxbs[1]\n# output\n4×4 Matrix{BigInt}:\n 5228766587057444  7771840027009430  7761130250550079  8285007185180178\n  725166959686948  3268240399638934  3257530623179583  3781407557809682\n 2976966773372196  5520040213324182  5509330436864831  6033207371494930\n 7480566400742692  1016440585953686  1005730809494335  1529607744124434","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"xbs[2]\n# output\n4×4 Matrix{BigInt}:\n 6432390138047189  6328265944608255   117225724257639  8118420346456289\n 1928790510676693  1824666317237759  4620825351628135  3614820719085793\n 4180590324361941  8580065758293503  2369025537942887  5866620532771041\n 8684189951732437  4076466130923007  6872625165313383  1363020905400545","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"BinaryToFloat64(xbs,rds_multiple)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Getting binary points with linear ordering is also supported. ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Reset!(ds) # resets rds_single and rds_multiple as well\nFirstLinearBinary(ds,2)\n# output\n4×4 Matrix{BigInt}:\n          0           0           0           0\n 2147483648  2147483648  2147483648  2147483648\n 1073741824  3221225472  3221225472  3221225472\n 3221225472  1073741824  1073741824  1073741824","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"FirstLinearBinary(rds_single,2)\n# output\n4×4 Matrix{BigInt}:\n 5228766587057444  6432390138047189  7771840027009430  6328265944608255\n  725166959686948  1928790510676693  3268240399638934  1824666317237759\n 7480566400742692  4180590324361941  1016440585953686  4076466130923007\n 2976966773372196  8684189951732437  5520040213324182  8580065758293503","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"xbs = FirstRLinearBinary(rds_multiple,2)\nxbs[1]\n# output\n4×4 Matrix{BigInt}:\n 5228766587057444  7771840027009430  7761130250550079  8285007185180178\n  725166959686948  3268240399638934  3257530623179583  3781407557809682\n 7480566400742692  1016440585953686  1005730809494335  1529607744124434\n 2976966773372196  5520040213324182  5509330436864831  6033207371494930","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"xbs[2]\n# output\n4×4 Matrix{BigInt}:\n 6432390138047189  6328265944608255   117225724257639  8118420346456289\n 1928790510676693  1824666317237759  4620825351628135  3614820719085793\n 8684189951732437  4076466130923007  6872625165313383  1363020905400545\n 4180590324361941  8580065758293503  2369025537942887  5866620532771041","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"These may be converted to floats as before. ","category":"page"},{"location":"tutorial/#IID-Standard-Uniform-Generator","page":"Tutorial","title":"IID Standard Uniform Generator","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"We provide an IID mathcalU01^s generator with the same API as Lattice and digital sequences. This is a wrapper around Random.MersenneTwister.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"For reproducibility, you may provide a seed.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"iiduseq = IIDU01Seq(3,7)\nNext(iiduseq,4)\n# output\n4×3 Matrix{Float64}:\n 0.812439   0.654977  0.489613\n 0.0787595  0.581591  0.258053\n 0.196465   0.193925  0.842951\n 0.66193    0.401352  0.635904","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Reset!(iiduseq)\nNext(iiduseq,4)\n# output\n4×3 Matrix{Float64}:\n 0.812439   0.654977  0.489613\n 0.0787595  0.581591  0.258053\n 0.196465   0.193925  0.842951\n 0.66193    0.401352  0.635904","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The seed gets used to construct a MersenneTwister. After","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using Random: MersenneTwister","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"You may pass also pass a MersenneTwister instance directly.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"iiduseq = IIDU01Seq(3,MersenneTwister(7))\nNext(iiduseq,4)\n# output\n4×3 Matrix{Float64}:\n 0.812439   0.654977  0.489613\n 0.0787595  0.581591  0.258053\n 0.196465   0.193925  0.842951\n 0.66193    0.401352  0.635904","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Providing neither a seed nor MersenneTwister uses MersenneTwister().","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"iiduseq = IIDU01Seq(3)\nsize(Next(iiduseq,4))\n# output\n(4, 3)","category":"page"},{"location":"tutorial/#Plotting","page":"Tutorial","title":"Plotting","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To save figures we need to ensure we are","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using CairoMakie","category":"page"},{"location":"tutorial/#Single-Projection","page":"Tutorial","title":"Single Projection","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"PLOTDIR = joinpath(@__DIR__,\"src/assets\")\nn = 2^6\nds = DigitalSeqB2G(3)\nfig = qmcscatter!(ds,n)\nsave(joinpath(PLOTDIR,\"basic.svg\"),fig)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#Extensibility","page":"Tutorial","title":"Extensibility","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"nvec = [1,2^6,2^7,2^8]\nfig = qmcscatter!(ds,nvec)\nsave(joinpath(PLOTDIR,\"extensibility.svg\"),fig)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#Multiple-Projections","page":"Tutorial","title":"Multiple Projections","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"dvec = [1 2; 1 3; 2 3]\nfig = qmcscatter!(ds,nvec,dvec)\nsave(joinpath(PLOTDIR,\"projections.svg\"),fig)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#Multiple-Randomizations","page":"Tutorial","title":"Multiple Randomizations","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"rds = RandomDigitalShift(DigitalSeqB2G(3),3)\nfig = qmcscatter!(rds,nvec,dvec)\nsave(joinpath(PLOTDIR,\"randomizations.svg\"),fig)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#Comparison-of-Sequences","page":"Tutorial","title":"Comparison of Sequences","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"iid = IIDU01Seq(3)\nrds = RandomDigitalShift(DigitalSeqB2G(3))\nrls = RandomShift(LatticeSeqB2(3))\nfig = qmcscatter!([1,2^6,2^7,2^8],[1 2],iid=iid,rds=rds,rls=rls)\nsave(joinpath(PLOTDIR,\"seq_comparison.svg\"),fig)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#MC-vs-QMC","page":"Tutorial","title":"MC vs QMC","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"We also need to be ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using Distributions\nusing LinearAlgebra","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"m = 16\nr = 100\nseed = 7\nn = 2^m\ns,mu = 7,-11.05684907978818\nrseqs = [IIDU01Seq(s,seed),RandomShift(LatticeSeqB2(s),r,seed),RandomDigitalShift(DigitalSeqB2G(s),r,seed)]\nxsets = [\n    [Next(rseqs[1],n) for k=1:r],\n    NextR(rseqs[2],n),\n    NextR(rseqs[3],n)]\nf(x::Vector{Float64}) = π^(s/2)*cos(norm(quantile.(Normal(),x)/sqrt(2)));\nf(x::Matrix{Float64}) = map(i->f(x[i,:]),1:size(x,1))\nfig = Figure(resolution=(800,500))\nax = Axis(fig[2,1],\n    xlabel = L\"$n$\",\n    ylabel = L\"$| \\hat{\\mu} - \\mu |$\",\n    yscale = log10,\n    xscale = log2)\nxlims!(ax,[1,n])\nfor k=1:size(xsets,1)\n    name,xs = rseqs[k].name,xsets[k]\n    ys = vcat(map(i->f(xs[i]),1:r)'...)\n    muhats = cumsum(ys,dims=2); for i=1:r muhats[i,:] = muhats[i,:]./[i for i=1:n] end \n    err = abs.(muhats.-mu)\n    pows2 = 2 .^ (0:m)\n    qlowerr = map(p2->quantile(err[:,p2],.35),pows2)\n    qmid = map(p2->quantile(err[:,p2],.5),pows2)\n    qhigherr = map(p2->quantile(err[:,p2],.65),pows2)\n    lines!(ax,pows2,qmid,color=JULIA4LOGOCOLORS[k],label=name,linewidth=3)\n    band!(pows2,qhigherr,qlowerr,color=(JULIA4LOGOCOLORS[k],.3))\nend\nfig[1,1] = Legend(fig,ax,framevisible=false,orientation=:horizontal)\nhidespines!(ax, :t, :r)\nsave(joinpath(PLOTDIR,\"mc_vs_qmc.svg\"),fig)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#Logo","page":"Tutorial","title":"Logo","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"nvec = [1,4,16,64]\nrds = RandomOwenScramble(DigitalSeqB2G(2),1,17)\nx = Next(rds,maximum(nvec))\nfig = Figure(resolution=(500,500),backgroundcolor=:transparent)\nax = Axis(fig[1,1],aspect=1,xticklabelsvisible=false,yticklabelsvisible=false,backgroundcolor=:transparent)\nqmcscatter!(ax,x,nvec)\nlimits!(ax,[-0.01,1.01],[-0.01,1.01])\nfor i=1:7 vlines!(ax,i/8,color=JULIA4LOGOCOLORS[3],alpha=1); hlines!(ax,i/8,color=JULIA4LOGOCOLORS[3],alpha=1) end \nfor i=1:3 vlines!(ax,i/4,color=JULIA4LOGOCOLORS[2],alpha=1); hlines!(ax,i/4,color=JULIA4LOGOCOLORS[2],alpha=1) end \nfor i=1:1 vlines!(ax,i/2,color=JULIA4LOGOCOLORS[1],alpha=1); hlines!(ax,i/2,color=JULIA4LOGOCOLORS[1],alpha=1) end \nhidespines!(ax); hidedecorations!(ax); hidexdecorations!(ax,grid = false); hideydecorations!(ax, ticks = false)\nsave(joinpath(PLOTDIR,\"logo.svg\"),fig)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"}]
}