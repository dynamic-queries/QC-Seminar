include("HHL.jl")
using Plots;

# See ReadMe for usage.
# Increase 1st parameter if you have a powerful machine
x_hhl,x_julia = driver(5,1);
plot(x_hhl,label=:"HHL");
plot!(x_julia,label=:"Julia");
