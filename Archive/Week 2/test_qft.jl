# Tests for routines in QFT.jl

include("qft.jl")

function test_qft_helpers()
    Ha = H();
    I_ = I();

    # >> Test 1 : H() Function
    k = 1/sqrt(2);
    @assert(Ha[1,1] == k);
    @assert(Ha[1,2] == k);
    @assert(Ha[2,1] == k);
    @assert(Ha[2,2] == -k);

    # >> Test 2 : I() Function

    @assert(I_[1,1] == 1.0);
    @assert(I_[1,2] == 0.0);
    @assert(I_[2,1] == 0.0);
    @assert(I_[2,2] == 1.0);

    # >> Test 3 : swap_function (scalar)
    a = 2;
    b = 3;
    a,b = swap(a,b);
    @assert(a==3);
    @assert(b==2);

    # >> Test 4 : flip vectors

    # Even case
    n = 6;
    y = rand(n);
    y_copy = copy(y);
    cswap(y);


    @assert(y[1] == y_copy[6]);
    @assert(y[2] == y_copy[5]);
    @assert(y[3] == y_copy[4]);

    # Odd case
    n = 5;
    y = rand(n);
    y_copy = copy(y);
    cswap(y);
    @assert(y[1] == y_copy[5]);
    @assert(y[2] == y_copy[4]);
    @assert(y[3] == y_copy[3]);
    println("All tests passed for the helper functions. \n");
end

function test_multiple_hadamard()
    Ha = H();
    Id = I();
    K = kron(Id,Ha,Id);
    # display(K);

    # println("\nTest candidate\n");
    B = HxN(3,2);
    # display(B);

    @assert(B==K);
    println("Multiple Hadamrd works like a charm!\n");
end


function test_controlled_rk()
    #=
        We shall compute a bunch for 3 qubit case
            q1  q2  q3
            c   t   i
            t   c   i
            i   c   t
            i   t   c
            c   i   t
            t   i   c
    =#
    n = 3;
    k =2;
    S = [1 0;0 exp(2*pi/2^(k))];
    Id = I();
    P0 = [1 0;0 0];
    P1 = [0 0;0 1];

    # For the lack of time, we only test the cases that we need.
    # This routine is not proven to be robust,
    # Functional for our implementation though.

    #Case 2
    O2 = kron(kron(Id,P0)+kron(S,P1),Id);
    A2 = Rx(n,k,2,1);
    @assert(O2==A2);

    # Case 6
    O3 = kron(kron(Id,Id),P0)+kron(kron(S,Id),P1);
    A3 = Rx(n,k,3,1);
    @assert(O3==A3);
    println("The Rkx function works adequately.");
end
