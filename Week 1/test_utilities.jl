include("Utilities.jl")
#= 
    1) test_discretize() :
    2) test_initial_state():  
=# 

function test_discretize() 
 
    #= **********Test case***************
    x0 = 0;
    x1 = 1;
    N = 11; 
    **Output :** 
    h = 1/10 = 0.1;
    x = 0,0,1.0,2.0,3....................1.0;
     ************************************ =#  
    x0 = 0;
    x1 = 1;
    N = 11; 

    x = discretize(x0,x1,N);
    # Length test 
    @assert(length(x) == N);
    # Test boundary conditions 
    @assert(x[1]==0);
    @assert(x[N]==1);
    # Test increment 
    @assert(x[2]-x[1] == 0.1); 
    println("All tests passed\n");
end

function test_initial_state() 
    #= **********Test case***************
    f(x) = 1; 
    x = dicretize(0,1,5);
    h = x[2]-x[1];
    tau = (h^2/2)/2;  // Cournant's condition 
    **Output**
    norm = 0.935414 
    T0 = {0.25,0.5,0.75};
    T_0_norm = {0.267261, 0.534522, 0.801784} 
    n_q = 4 
    ************************************=#
    x = discretize(0,1,5);
    h = x[2] - x[1];
    f(x) = x; 
    tau = 0; 
    n_q,T0 = initial_state(x,f,tau);
    # Assertions 
    # Functional Values 
    @assert(isapprox(T0[1],0.267261,atol=1e-1));
    @assert(isapprox(T0[2],0.534522,atol=1e-1));
    @assert(isapprox(T0[3],0.801784,atol=1e-1));
    println("All tests passed\n"); 
end
