```@meta
DocTestSetup = quote
    using QMCGenerators
    using CairoMakie
    using Random: MersenneTwister
    using Distributions
    using LinearAlgebra
end
```

# Tutorial

```@contents
Pages = ["tutorial.md"]
Depth = 4
```

To begin, install this package with 

```julia 
Pkg.add("QMCGenerators")
```

and then import via

```julia
using QMCGenerators
```

## Common Usage

Often we just want to generate a single randomized quasi-random sequence. For example, to generate the first 4 points of a 3 dimensional digital net with 1 random digital shift with seed 7

```jldoctest
Next(RandomDigitalShift(DigitalSeqB2G(3),1,7),4)
# output
4×3 Matrix{Float64}:
 0.648627  0.0158008  0.971646
 0.148627  0.515801   0.471646
 0.398627  0.265801   0.721646
 0.898627  0.765801   0.221646
```

A similar API is available for randomly shifted Lattices:

```jldoctest
Next(RandomShift(LatticeSeqB2(3),1,7),4)
# output
4×3 Matrix{Float64}:
 0.812439   0.654977  0.489613
 0.312439   0.154977  0.989613
 0.0624386  0.404977  0.239613
 0.562439   0.904977  0.739613
```

!!! warning
    While not strictly enforced, sample sizes should be powers of two to achieve full coverage of $[0,1]^s$

In the following sections we always supply a seed for reproducibility. Supplying a seed requires you also supply the number of randomizations as done above. However, if you do not wish to seed, you can simply supply the number of randomizations. 

```jldoctest tut_no_seed
rls = RandomDigitalShift(DigitalSeqB2G(12),2)
xs = NextR(rls,2^7)
size(xs)
# output
(2,)
```
```jldoctest tut_no_seed
size(xs[1]) == size(xs[2]) == (2^7,12)
# output
true
```

Moreover, you may use only one randomization without a seed with the simplified API
```jldoctest
rls = RandomDigitalShift(DigitalSeqB2G(52))
x = Next(rls,2^14)
size(x)
# output
(16384, 52)
```

The same API simplifications holds for Lattices. 

## Structure and Functions

### Unrandomized Sequences

!!! warning
    It is highly recommended you randomize sequences. The first point of an unrandomized sequence is $0 \in [0,1]^s$ which will be transformed to an infinite values in many functions e.g. those composed with taking the inverse CDF of a $\mathcal{N}(0,1)$ as done in the section on [Quasi-Monte Carlo](@ref). 


Let's start by defining a 5 dimensional (unrandomized) digital sequence and generating the first 4 points.

```jldoctest tut_ds
ds = DigitalSeqB2G(5)
Next(ds,4)
# output
4×5 Matrix{Float64}:
 0.0   0.0   0.0   0.0   0.0
 0.5   0.5   0.5   0.5   0.5
 0.75  0.25  0.25  0.25  0.75
 0.25  0.75  0.75  0.75  0.25
```

We can then generate the next 4 points in the sequence with 

```jldoctest tut_ds
Next(ds,4)
# output
4×5 Matrix{Float64}:
 0.375  0.375  0.625  0.875  0.375
 0.875  0.875  0.125  0.375  0.875
 0.625  0.125  0.875  0.625  0.625
 0.125  0.625  0.375  0.125  0.125
```

To reset the generator use

```jldoctest tut_ds
Reset!(ds)
Next(ds,4)
# output
4×5 Matrix{Float64}:
 0.0   0.0   0.0   0.0   0.0
 0.5   0.5   0.5   0.5   0.5
 0.75  0.25  0.25  0.25  0.75
 0.25  0.75  0.75  0.75  0.25
```

Let's reset once more before continuing 

```jldoctest tut_ds; output = false
Reset!(ds)
# output
```

These functions can also be applied to Lattices 

```jldoctest tut_ls
ls = LatticeSeqB2(5)
Next(ls,4)
# output
4×5 Matrix{Float64}:
 0.0   0.0   0.0   0.0   0.0
 0.5   0.5   0.5   0.5   0.5
 0.25  0.75  0.75  0.25  0.25
 0.75  0.25  0.25  0.75  0.75
```
```jldoctest tut_ls
Next(ls,4)
# output
4×5 Matrix{Float64}:
 0.125  0.375  0.375  0.125  0.125
 0.625  0.875  0.875  0.625  0.625
 0.375  0.125  0.125  0.375  0.375
 0.875  0.625  0.625  0.875  0.875
```
```jldoctest tut_ls; output = false
Reset!(ls)
# output
```

###  Independent Sequence Randomizations

Getting a single randomization was shown in the [Common Usage](@ref) section. Support is also available for multiple independent randomizations. For instance, we can generating 2 independent randomizations with seed 11 and get the next 4 points in each sequence via  

```jldoctest tut_ds
rds = RandomDigitalShift(ds,2,11)
xs = NextR(rds,4)
typeof(xs)
# output
Vector{Matrix{Float64}} (alias for Array{Array{Float64, 2}, 1})
```
```jldoctest tut_ds
xs[1]
# output
4×5 Matrix{Float64}:
 0.50262     0.64885  0.771298   0.211364  0.206083
 0.00261977  0.14885  0.271298   0.711364  0.706083
 0.25262     0.89885  0.521298   0.461364  0.956083
 0.75262     0.39885  0.0212976  0.961364  0.456083
```
```jldoctest tut_ds
xs[2]
# output
4×5 Matrix{Float64}:
 0.583259   0.861213  0.781133  0.458726  0.93243
 0.0832589  0.361213  0.281133  0.958726  0.43243
 0.333259   0.611213  0.531133  0.208726  0.18243
 0.833259   0.111213  0.031133  0.708726  0.68243
```

As with unrandomized sequences, we can get the next 4 points and then reset the generator with 

```jldoctest tut_ds
xs = NextR(rds,4)
xs[1]
# output
4×5 Matrix{Float64}:
 0.87762  0.77385  0.396298  0.836364   0.331083
 0.37762  0.27385  0.896298  0.336364   0.831083
 0.12762  0.52385  0.146298  0.586364   0.581083
 0.62762  0.02385  0.646298  0.0863639  0.0810831
```
```jldoctest tut_ds
xs[2]
# output
4×5 Matrix{Float64}:
 0.958259  0.736213  0.406133  0.583726   0.55743
 0.458259  0.236213  0.906133  0.0837258  0.0574296
 0.208259  0.986213  0.156133  0.833726   0.30743
 0.708259  0.486213  0.656133  0.333726   0.80743
```
```jldoctest tut_ds; output = false
Reset!(rds)
# output
```

Similarly for Lattices 

```jldoctest tut_ls
rls = RandomShift(ls,2,11)
xr = NextR(rls,4)
xr[1]
# output
4×5 Matrix{Float64}:
 0.408169  0.860346  0.338226   0.58282    0.961567
 0.908169  0.360346  0.838226   0.0828199  0.461567
 0.658169  0.610346  0.0882258  0.83282    0.211567
 0.158169  0.110346  0.588226   0.33282    0.711567
```
```jldoctest tut_ls
xr[2]
# output
4×5 Matrix{Float64}:
 0.515953   0.634856  0.376913  0.695549  0.187747
 0.0159529  0.134856  0.876913  0.195549  0.687747
 0.765953   0.384856  0.126913  0.945549  0.437747
 0.265953   0.884856  0.626913  0.445549  0.937747
```
```jldoctest tut_ls; output = false
Reset!(rls)
# output
```

## Advanced Features

### Alternative Generating Matrices and Vectors

#### Pregenerated

We include many alternative generating matrices for digital sequences in [this directory](https://bitbucket.org/dnuyens/qmc-generators/src/master/DIGSEQ/). To use an alternative, simply supply the relative path 

```jldoctest 
ds = DigitalSeqB2G(3,"sobolmats/sobol_alpha2_Bs64.col")
Next(ds,4)
# output
4×3 Matrix{Float64}:
 0.0     0.0     0.0
 0.75    0.75    0.75
 0.6875  0.1875  0.9375
 0.4375  0.9375  0.1875
```

Alternative Lattice generating vectors are available in [this directory](https://bitbucket.org/dnuyens/qmc-generators/src/master/LATSEQ/). For Lattices, after supplying the path you also need to pass the $m$ value in the file name

```jldoctest
ls = LatticeSeqB2(3,"exod8_base2_m13.txt",13)
Next(ls,4)
# output
4×3 Matrix{Float64}:
 0.0   0.0   0.0
 0.5   0.5   0.5
 0.25  0.25  0.75
 0.75  0.75  0.25
```

#### User Defined

One may supply their own generating matrix to construct a base 2 digital sequence, for example

```jldoctest tut_ds_custom_matrix
m = 5
C1 = [BigInt(2^i) for i=0:(m-1)]
C2 = [BigInt(1) for i=1:m]
for i in 2:m C2[i] = (C2[i-1] << 1) ⊻ C2[i-1] end
generating_matrix = vcat(C1',C2')
# output
2×5 Matrix{BigInt}:
 1  2  4   8  16
 1  3  5  15  17
```
```jldoctest tut_ds_custom_matrix
ds = DigitalSeqB2G(2,generating_matrix)
Next(ds,4)
# output
4×2 Matrix{Float64}:
 0.0   0.0
 0.5   0.5
 0.75  0.25
 0.25  0.75
```

For base 2 Lattices, you may supply the generating vector followed by $m$ where $2^m$ is the maximum number of supported points 

```jldoctest
generating_vector = BigInt[1,433461,315689]
m = 20
ls = LatticeSeqB2(3,generating_vector,m)
Next(ls,4)
# output
4×3 Matrix{Float64}:
 0.0   0.0   0.0
 0.5   0.5   0.5
 0.25  0.25  0.25
 0.75  0.75  0.75
```

### Linear Ordering

By default, digital sequences are generated in Gray code order. One may generate the first $2^m$ points in linear order via 

```jldoctest tut_ds_order
m = 3
n = 2^m
ds = DigitalSeqB2G(4)
FirstLinear(ds,m)
# output
8×4 Matrix{Float64}:
 0.0    0.0    0.0    0.0
 0.5    0.5    0.5    0.5
 0.25   0.75   0.75   0.75
 0.75   0.25   0.25   0.25
 0.125  0.625  0.375  0.125
 0.625  0.125  0.875  0.625
 0.375  0.375  0.625  0.875
 0.875  0.875  0.125  0.375
```

Compare to the original ordering 

```jldoctest tut_ds_order
Next(ds,n)
# output
8×4 Matrix{Float64}:
 0.0    0.0    0.0    0.0
 0.5    0.5    0.5    0.5
 0.75   0.25   0.25   0.25
 0.25   0.75   0.75   0.75
 0.375  0.375  0.625  0.875
 0.875  0.875  0.125  0.375
 0.625  0.125  0.875  0.625
 0.125  0.625  0.375  0.125
```

Similarly, Lattices are by default generated in extensible ordering. A linear ordering is also available 

```jldoctest tut_ls_order
m = 3
n = 2^m 
ls = LatticeSeqB2(4)
FirstLinear(ls,m)
# output
8×4 Matrix{Float64}:
 0.0    0.0    0.0    0.0
 0.125  0.375  0.375  0.125
 0.25   0.75   0.75   0.25
 0.375  0.125  0.125  0.375
 0.5    0.5    0.5    0.5
 0.625  0.875  0.875  0.625
 0.75   0.25   0.25   0.75
 0.875  0.625  0.625  0.875
```

Compare to the original ordering 

```jldoctest tut_ls_order
Next(ls,n)
# output
8×4 Matrix{Float64}:
 0.0    0.0    0.0    0.0
 0.5    0.5    0.5    0.5
 0.25   0.75   0.75   0.25
 0.75   0.25   0.25   0.75
 0.125  0.375  0.375  0.125
 0.625  0.875  0.875  0.625
 0.375  0.125  0.125  0.375
 0.875  0.625  0.625  0.875
```

Linear order for randomized sequences has expected syntax 

```jldoctest tut_ds_order_rand
ds = DigitalSeqB2G(4)
rds = RandomDigitalShift(ds,1,17)
FirstLinear(rds,2)
# output
4×4 Matrix{Float64}:
 0.971405  0.428212  0.121174  0.444364
 0.471405  0.928212  0.621174  0.944364
 0.721405  0.678212  0.871174  0.694364
 0.221405  0.178212  0.371174  0.194364
```
```jldoctest tut_ds_order_rand
rds = RandomDigitalShift(ds,2,17)
xs = FirstRLinear(rds,2)
xs[1]
# output
4×4 Matrix{Float64}:
 0.971405  0.121174  0.0747568  0.82345
 0.471405  0.621174  0.574757   0.32345
 0.721405  0.871174  0.824757   0.07345
 0.221405  0.371174  0.324757   0.57345
```

```jldoctest tut_ds_order_rand
xs[2]
# output
4×4 Matrix{Float64}:
 0.428212  0.444364  0.320031   0.705377
 0.928212  0.944364  0.820031   0.205377
 0.178212  0.694364  0.570031   0.455377
 0.678212  0.194364  0.0700311  0.955377
```

The same functions are available for randomly shifted lattices. 

### Binary Functions for Digital Sequences

For digital sequences, we sometimes want the binary representation of points. We can get the binary representations as integers and then convert them to their floating point values as follows 

```jldoctest tut_ds_binary
ds = DigitalSeqB2G(4)
xb = NextBinary(ds,4)
# output
4×4 Matrix{BigInt}:
          0           0           0           0
 2147483648  2147483648  2147483648  2147483648
 3221225472  1073741824  1073741824  1073741824
 1073741824  3221225472  3221225472  3221225472
```
```jldoctest tut_ds_binary
BinaryToFloat64(xb,ds)
# output
4×4 Matrix{Float64}:
 0.0   0.0   0.0   0.0
 0.5   0.5   0.5   0.5
 0.75  0.25  0.25  0.25
 0.25  0.75  0.75  0.75
```

This is also compatible with randomized digital sequences 

```jldoctest tut_ds_binary
Reset!(ds)
rds_single = RandomDigitalShift(ds,1,11)
xb = NextBinary(rds_single,4)
# output
4×4 Matrix{BigInt}:
 4527196440290335  5253528859752444  5844321400240638  7757115741280515
   23596812919839   749929232381948  1340721772870142  3253516113910019
 2275396626605087  7505328673437692  8096121213925886  5505315927595267
 6778996253975583  3001729046067196  3592521586555390  1001716300224771
```
```jldoctest tut_ds_binary
BinaryToFloat64(xb,rds_single)
# output
4×4 Matrix{Float64}:
 0.50262     0.583259   0.64885  0.861213
 0.00261977  0.0832589  0.14885  0.361213
 0.25262     0.833259   0.89885  0.611213
 0.75262     0.333259   0.39885  0.111213
```
```jldoctest tut_ds_binary
Reset!(ds)
rds_multiple = RandomDigitalShift(ds,2,11)
xbs = NextRBinary(rds_multiple,4)
xbs[1]
# output
4×4 Matrix{BigInt}:
 4527196440290335  5844321400240638  6947231307756753  1903796900808484
   23596812919839  1340721772870142  2443631680386257  6407396528178980
 2275396626605087  8096121213925886  4695431494071505  4155596714493732
 6778996253975583  3592521586555390   191831866701009  8659196341864228
```
```jldoctest tut_ds_binary
xbs[2]
# output
4×4 Matrix{BigInt}:
 5253528859752444  7757115741280515  7035820915317184  4131834350673793
  749929232381948  3253516113910019  2532221287946688  8635433978044289
 3001729046067196  5505315927595267  4784021101631936  1880034536988545
 7505328673437692  1001716300224771   280421474261440  6383634164359041
```
```jldoctest tut_ds_binary; output = false
BinaryToFloat64(xbs,rds_multiple)
# output
2-element Vector{Matrix{Float64}}:
 [0.5026197725011378 0.6488500181856691 0.7712976155268284 0.21136391534875942; 0.0026197725011377626 0.14885001818566912 0.2712976155268284 0.7113639153487594; 0.25261977250113776 0.8988500181856691 0.5212976155268284 0.4613639153487594; 0.7526197725011378 0.3988500181856691 0.021297615526828406 0.9613639153487594]
 [0.5832588700629908 0.8612128500651867 0.7811330377324381 0.4587257630055178; 0.08325887006299082 0.36121285006518666 0.28113303773243814 0.9587257630055178; 0.3332588700629908 0.6112128500651867 0.5311330377324381 0.20872576300551782; 0.8332588700629908 0.11121285006518666 0.03113303773243814 0.7087257630055178]
```

Getting binary points with linear ordering is also supported. 

```jldoctest tut_ds_binary
Reset!(ds) # resets rds_single and rds_multiple as well
FirstLinearBinary(ds,2)
# output
4×4 Matrix{BigInt}:
          0           0           0           0
 2147483648  2147483648  2147483648  2147483648
 1073741824  3221225472  3221225472  3221225472
 3221225472  1073741824  1073741824  1073741824
```
```jldoctest tut_ds_binary
FirstLinearBinary(rds_single,2)
# output
4×4 Matrix{BigInt}:
 4527196440290335  5253528859752444  5844321400240638  7757115741280515
   23596812919839   749929232381948  1340721772870142  3253516113910019
 6778996253975583  3001729046067196  3592521586555390  1001716300224771
 2275396626605087  7505328673437692  8096121213925886  5505315927595267
```
```jldoctest tut_ds_binary
xbs = FirstRLinearBinary(rds_multiple,2)
xbs[1]
# output
4×4 Matrix{BigInt}:
 4527196440290335  5844321400240638  6947231307756753  1903796900808484
   23596812919839  1340721772870142  2443631680386257  6407396528178980
 6778996253975583  3592521586555390   191831866701009  8659196341864228
 2275396626605087  8096121213925886  4695431494071505  4155596714493732
```
```jldoctest tut_ds_binary
xbs[2]
# output
4×4 Matrix{BigInt}:
 5253528859752444  7757115741280515  7035820915317184  4131834350673793
  749929232381948  3253516113910019  2532221287946688  8635433978044289
 7505328673437692  1001716300224771   280421474261440  6383634164359041
 3001729046067196  5505315927595267  4784021101631936  1880034536988545
```

These may be converted to floats as before. 


### IID Standard Uniform Generator

We provide an IID $\mathcal{U}[0,1]^s$ generator with the same API as Lattice and digital sequences. This is a wrapper around [`Random.MersenneTwister`](https://docs.julialang.org/en/v1/stdlib/Random/#Random.MersenneTwister).

For reproducibility, you may provide a seed.

```jldoctest iidu
iiduseq = IIDU01Seq(3,7)
Next(iiduseq,4)
# output
4×3 Matrix{Float64}:
 0.812439   0.654977  0.489613
 0.0787595  0.581591  0.258053
 0.196465   0.193925  0.842951
 0.66193    0.401352  0.635904
```
```jldoctest iidu
Reset!(iiduseq)
Next(iiduseq,4)
# output
4×3 Matrix{Float64}:
 0.812439   0.654977  0.489613
 0.0787595  0.581591  0.258053
 0.196465   0.193925  0.842951
 0.66193    0.401352  0.635904
```

The seed gets used to construct a `MersenneTwister`. After

```julia
using Random: MersenneTwister
```

You may pass also pass a `MersenneTwister` instance directly.

```jldoctest iidu
iiduseq = IIDU01Seq(3,MersenneTwister(7))
Next(iiduseq,4)
# output
4×3 Matrix{Float64}:
 0.812439   0.654977  0.489613
 0.0787595  0.581591  0.258053
 0.196465   0.193925  0.842951
 0.66193    0.401352  0.635904
```

Providing neither a seed nor `MersenneTwister` uses `MersenneTwister()`.

```jldoctest iidu
iiduseq = IIDU01Seq(3)
size(Next(iiduseq,4))
# output
(4, 3)
```

## Plotting

To save figures we need to ensure we are

```julia 
using CairoMakie
```

### Single Projection

```jldoctest plots; output = false
PLOTDIR = joinpath(@__DIR__,"src/assets")
n = 2^6
ds = DigitalSeqB2G(3)
fig = qmcscatter!(ds,n)
save(joinpath(PLOTDIR,"basic.svg"),fig)
# output
CairoMakie.Screen{SVG}
```

![image](./assets/basic.svg)

### Extensibility 

```jldoctest plots; output = false
nvec = [1,2^6,2^7,2^8]
fig = qmcscatter!(ds,nvec)
save(joinpath(PLOTDIR,"extensibility.svg"),fig)
# output
CairoMakie.Screen{SVG}
```

![image](./assets/extensibility.svg)

### Multiple Projections

```jldoctest plots; output = false
dvec = [1 2; 1 3; 2 3]
fig = qmcscatter!(ds,nvec,dvec)
save(joinpath(PLOTDIR,"projections.svg"),fig)
# output
CairoMakie.Screen{SVG}
```

![image](./assets/projections.svg)

### Multiple Randomizations

```jldoctest plots; output = false
rds = RandomDigitalShift(DigitalSeqB2G(3),3)
fig = qmcscatter!(rds,nvec,dvec)
save(joinpath(PLOTDIR,"randomizations.svg"),fig)
# output
CairoMakie.Screen{SVG}
```

![image](./assets/randomizations.svg)

### Comparison of Sequences

```jldoctest plots; output = false
iid = IIDU01Seq(3)
rds = RandomDigitalShift(DigitalSeqB2G(3))
rls = RandomShift(LatticeSeqB2(3))
fig = qmcscatter!([1,2^6,2^7,2^8],[1 2],iid=iid,rds=rds,rls=rls)
save(joinpath(PLOTDIR,"seq_comparison.svg"),fig)
# output
CairoMakie.Screen{SVG}
```

![image](./assets/seq_comparison.svg)

### MC vs QMC

We also need to be 

```julia 
using Distributions
using LinearAlgebra
```

```jldoctest plots; output = false
m = 16
r = 100
seed = 7
n = 2^m
s,mu = 7,-11.05684907978818
rseqs = [IIDU01Seq(s,seed),RandomShift(LatticeSeqB2(s),r,seed),RandomDigitalShift(DigitalSeqB2G(s),r,seed)]
xsets = [
    [Next(rseqs[1],n) for k=1:r],
    NextR(rseqs[2],n),
    NextR(rseqs[3],n)]
f(x::Vector{Float64}) = π^(s/2)*cos(norm(quantile.(Normal(),x)/sqrt(2)));
f(x::Matrix{Float64}) = map(i->f(x[i,:]),1:size(x,1))
fig = Figure(resolution=(800,500))
ax = Axis(fig[2,1],
    xlabel = L"$n$",
    ylabel = L"$| \hat{\mu} - \mu |$",
    yscale = log10,
    xscale = log2)
xlims!(ax,[1,n])
for k=1:size(xsets,1)
    name,xs = rseqs[k].name,xsets[k]
    ys = vcat(map(i->f(xs[i]),1:r)'...)
    muhats = cumsum(ys,dims=2); for i=1:r muhats[i,:] = muhats[i,:]./[i for i=1:n] end 
    err = abs.(muhats.-mu)
    pows2 = 2 .^ (0:m)
    qlowerr = map(p2->quantile(err[:,p2],.35),pows2)
    qmid = map(p2->quantile(err[:,p2],.5),pows2)
    qhigherr = map(p2->quantile(err[:,p2],.65),pows2)
    lines!(ax,pows2,qmid,color=JULIA4LOGOCOLORS[k],label=name,linewidth=3)
    band!(pows2,qhigherr,qlowerr,color=(JULIA4LOGOCOLORS[k],.3))
end
fig[1,1] = Legend(fig,ax,framevisible=false,orientation=:horizontal)
hidespines!(ax, :t, :r)
save(joinpath(PLOTDIR,"mc_vs_qmc.svg"),fig)
# output
CairoMakie.Screen{SVG}
```

![image](./assets/mc_vs_qmc.svg)

### Logo 

```jldoctest plots; output = false
nvec = [1,4,16,64]
rds = RandomDigitalShift(DigitalSeqB2G(2),1,17)
x = Next(rds,maximum(nvec))
fig = Figure(resolution=(500,500),backgroundcolor=:transparent)
ax = Axis(fig[1,1],aspect=1,xticklabelsvisible=false,yticklabelsvisible=false,backgroundcolor=:transparent)
qmcscatter!(ax,x,nvec)
limits!(ax,[-0.01,1.01],[-0.01,1.01])
for i=1:7 vlines!(ax,i/8,color=JULIA4LOGOCOLORS[3],alpha=1); hlines!(ax,i/8,color=JULIA4LOGOCOLORS[3],alpha=1) end 
for i=1:3 vlines!(ax,i/4,color=JULIA4LOGOCOLORS[2],alpha=1); hlines!(ax,i/4,color=JULIA4LOGOCOLORS[2],alpha=1) end 
for i=1:1 vlines!(ax,i/2,color=JULIA4LOGOCOLORS[1],alpha=1); hlines!(ax,i/2,color=JULIA4LOGOCOLORS[1],alpha=1) end 
hidespines!(ax); hidedecorations!(ax); hidexdecorations!(ax,grid = false); hideydecorations!(ax, ticks = false)
save(joinpath(PLOTDIR,"logo.svg"),fig)
# output
CairoMakie.Screen{SVG}
```

![image](./assets/logo.svg)
