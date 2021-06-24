#=
Task 2 : To compute the eigenvalues of the Simulation Matrix  using Phase Estimation

1) Implementations to Phase Estimation Routines

2) Inverse Fourier Transform and Forward Fourier Transform implementations

    i) Hadamard Gate
   ii) Conttrolled Rotation
  iii) Swap Gate

=#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function H()
    A = Array{Float64}(undef,2,2);
    A[1,1] = 1; A[1,2] = 1; A[2,1] = 1; A[2,2] = -1;
    return 1/sqrt(2) * A;
end
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function I()
    I = Array{Float64}(undef,2,2);
    I[1,1] = 1.0; I[2,2] = 1.0;I[1,2] = 0.0; I[2,1] =0.0;
    return I;
end
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function CU(n)
    return [1 0; 0 exp(2*pi*1im/2^n)]
end
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function swap(a,b)
    temp = a;
    a = b;
    b = temp;
    return a,b;
end
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function cswap(y)
    #=
    This is a classical swap operation.
    It could be done using quantum gates equivalently.
    However this could lead to needless computations.
    When swapping gates we have two cases :
    1) If n is even we swap the extreme gates but the one in the middle
    2) If n is odd we swap all the extreme gates
    =#
    n = length(y);
    for i=1:Int16(ceil(n/2))
            y[i],y[n-i+1] = swap(y[i],y[n-i+1]);
    end
 end
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function HxN(n,i)
    #=
    #---------------------------------------------------------------------------
    How to use this : n is the number of qubits. i is the location of H operator.
    #---------------------------------------------------------------------------
    Scope for the future : Expand functionality to all Operators, by including an argument.
    =#

    # Initializing the operator to the first matrix
    #---------------------------------------------------------------------------
    H_ = H();
    I_ = I();
    #---------------------------------------------------------------------------
    if(i==1)
        Hxn = H_;
    else
        Hxn = I_;
    end
    #---------------------------------------------------------------------------
    # Moving on
    for j=2:n
        if (j==i)
            Hxn = kron(Hxn,H_);
        else
            Hxn = kron(Hxn,I_);
        end
    end
    #---------------------------------------------------------------------------
    return Hxn
end

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function Rx(n,k,control,target)
    # Since we are dealing with IQFT, the matrix is hardcoded. However this could be done
    # for any matrix.
    #=
        kRot gate has a control gate and a target gate.
        As we know a controlled operation is done as follows.
        (P0 x I) + (P1 x S)
    =#
    @assert(control!=target);

    S = [1 0;0 exp(2*pi/2^(k))];
    Id = I();
    P0 = [1 0;0 0];
    P1 = [0 0;0 1];

    O1 = Id;
    O2 = Id;
    for i=1:n
        if (i==1)
            if (target == 1)
                O1 = Id;
                O2 = S;
            else
                O1 = P0;
                O2 = P1;
            end
            continue;
        end
        if (i==target)
            O1 = kron(O1,Id);
            O2 = kron(O2,S);
        elseif (i==control)
            O1 = kron(O1,P0);
            O2 = kron(O2,P1);
        else
            O1 = kron(O1,Id);
            O2 = kron(O2,Id);
        end
    end
    O = O1 + O2;
    return O;
end
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function IQFT(y)
    #=
    Let y be an arbitrary state vector
    If y is comprised of n qubits => length(y) = 2^n => n=log(2,length(y));
    =#
    n_q = log(2,length(y));
    y = cswap(y);
    for i=n_q:-1:1

    end
end
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
