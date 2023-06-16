```@meta
DocTestSetup = quote
    using QMCGenerators
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
5-element Vector{BigInt}:
 0
 0
 0
 0
 0
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
-1
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
5-element Vector{BigInt}:
 0
 0
 0
 0
 0
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
-1
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
 2158735485  2505077772  2786789608  3698881025
   11251837   357594124   639305960  1551397377
 1084993661  3578819596  3860531432  2625139201
 3232477309  1431335948  1713047784   477655553
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
 2158735485  2786789608  3312698034   907801103
   11251837   639305960  1165214386  3055284751
 1084993661  3860531432  2238956210  1981542927
 3232477309  1713047784    91472562  4129026575
```
```jldoctest tut_ds_binary
xbs[2]
# output
4×4 Matrix{BigInt}:
 2505077772  3698881025  3354940850  1970212149
  357594124  1551397377  1207457202  4117695797
 1431335948  2625139201  2281199026   896470325
 3578819596   477655553   133715378  3043953973
```
```jldoctest tut_ds_binary; output = false
BinaryToFloat64(xbs,rds_multiple)
# output
2-element Vector{Matrix{Float64}}:
 [0.5026197724509984 0.648850018158555 0.7712976154871285 0.2113639151211828; 0.0026197724509984255 0.14885001815855503 0.2712976154871285 0.7113639151211828; 0.2526197724509984 0.898850018158555 0.5212976154871285 0.4613639151211828; 0.7526197724509984 0.39885001815855503 0.021297615487128496 0.9613639151211828]
 [0.5832588700577617 0.861212849849835 0.7811330375261605 0.4587257627863437; 0.08325887005776167 0.36121284984983504 0.2811330375261605 0.9587257627863437; 0.33325887005776167 0.611212849849835 0.5311330375261605 0.2087257627863437; 0.8332588700577617 0.11121284984983504 0.03113303752616048 0.7087257627863437]
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
 2158735485  2505077772  2786789608  3698881025
   11251837   357594124   639305960  1551397377
 3232477309  1431335948  1713047784   477655553
 1084993661  3578819596  3860531432  2625139201
```
```jldoctest tut_ds_binary
xbs = FirstRLinearBinary(rds_multiple,2)
xbs[1]
# output
4×4 Matrix{BigInt}:
 2158735485  2786789608  3312698034   907801103
   11251837   639305960  1165214386  3055284751
 3232477309  1713047784    91472562  4129026575
 1084993661  3860531432  2238956210  1981542927
```
```jldoctest tut_ds_binary
xbs[2]
# output
4×4 Matrix{BigInt}:
 2505077772  3698881025  3354940850  1970212149
  357594124  1551397377  1207457202  4117695797
 3578819596   477655553   133715378  3043953973
 1431335948  2625139201  2281199026   896470325
```


These may be converted to floats as before. 




