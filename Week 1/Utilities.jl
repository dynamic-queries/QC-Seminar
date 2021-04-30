using LinearAlgebra 


#= 
    Task 1 : 
            Function to represent the initial normalized state vector in  Cn 
    
    1) 
    Function name : intial_state 
    Requirements  : 
                    1) Compute the number of qubits required to represent the domain
                    2) Calculate the intial functional approximation in the domain 
                    3) Compure norm of the functional set
                    4) Normalize the initial functional set 
    Params 1      : x,f(x)            // Domain for discretization and intial state function  
    Output        : n, T0             // Number of qubits required, Initial state in Cn 

    2) 
    Function name : discretize 
    Requirements  : Discretize the domain using the passed parameters 
    Params 1      : x_initial , x_final , N    // N := number of discretization points 
    Output        : x                          // Domain of discretization

    3) 
=# 

function discretize(x_intial,x_final,N)
    x = LinRange(x_intial,x_final,N);
    return x; 
end

function f(x)
    return x;
end

function initial_state(x,f,tau)
    n = length(x); 
    # But in the discretization equation Ax = b we do not account boundary values. So N ~= N-2;  
    n_qubits = log(2,n-1);
    h = x[2]-x[1]; 
    # Compute f(x_i) 
    T_0 = f(x[2:n-1]);  
    # Add boundary conditions on to b 
    temp1 = zeros(n-2);
    temp1[1] = tau/h^2 * f(x[1]);
    temp1[2] = tau/h^2 * f(x[n]);
    temp2 = temp1 + T_0; 
    T_0 = temp2;
    # Norm of T_0 
    T_0 = T_0 / norm(T_0); 
    return n_qubits,T_0; 
end

