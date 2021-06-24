### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ cf4de4d2-b32f-42b8-828e-aa9a15e0924c
begin
	using Yao
	using YaoExtensions
	using Plots
	using BitBasis
	using LinearAlgebra
	using PlutoUI
end

# ╔═╡ 3b41ca30-d458-11eb-1db5-151b58667f27
function compute_θ(nreg,scale,t,ntot)
	# We estimate angles corresponding to every element of the state vector
	circuit = chain(ntot);
	l =  2^nreg;
	angles = zeros(l-1);
	for i=1:l-1
		val = (scale*l)/i;
		if (val ≈ 1)
			angles[i] = 2*asin(val);
		elseif (val < 1)
			angles[i] = 2*asin(val);
		else
			angles[i] = 0 ;
		end
	# We now have the list of all possible angles that can arise
	# This means that there should be l possible gate implementations
	control_vec = Vector(2:(nreg+1));
	for j=1:(nreg)
		if(readbit(i,j)==0)
			control_vec[j] = -control_vec[j];
		end
	end
	push!(circuit,control(control_vec,1=>Ry(angles[i]/t)));
	end
	return circuit;
end

# ╔═╡ 5fe8135e-7f83-4b20-b36e-b593153e333c
function transform(nreg,scaling)
	ntot = nreg+1
	circuit = chain(ntot);
	for j = 2:ntot
		sub_circuit = compute_θ(nreg,scaling,j-1,ntot);
		circuit = chain(circuit,sub_circuit);
		nreg = nreg-1;
	end
	return circuit;
end

# ╔═╡ 5d8f93fa-2598-4ae1-8d67-9befe6385e7d
function PE(neig,nvec,U)
	ntotal = Int(neig + nvec);
	uniform_superposition =  repeat(ntotal,H,1:neig);
	quantum_simulation = chain(ntotal);
	U = GeneralMatrixBlock(U);
	for i=1:neig
		push!(quantum_simulation,control(ntotal,(i,),[neig+1:ntotal...,]=>U));
		if(i!=neig)
			U = matblock(mat(U)*mat(U));
		end
	end
	inverse = subroutine(ntotal,qft(neig)',[1:neig...,]);
	circuit = chain(uniform_superposition,quantum_simulation,inverse);
	return circuit
end

# ╔═╡ caff44b6-72b8-4eca-8e79-6b47f163612e
function isValid(A,U,b)
	@assert(A'≈ A);
	@assert(U*U' ≈ I);
	@assert(norm(b) ≈ 1);
end

# ╔═╡ 7de7d7af-66a6-477d-a43b-3b394bcdf3ac
begin
	struct CRot{N, NC, T} <: PrimitiveBlock{N}
	    cbits::Vector{Int}
	    ibit::Int
	    C_value::T
	    CRot{N}(cbits::Vector{Int}, ibit::Int, C_value::T) where {N, T} = new{N, length(cbits), T}(cbits, ibit, C_value)
	end

	@inline function rotmat(λ::Real, C_value::Real)
	    b = C_value/λ
	    a = sqrt(1-b^2)
	    a, -b, b, a
	end
	#=
	Reused under CC License from Yao.jl
	=#
	function YaoBlocks._apply!(reg::ArrayReg, hr::CRot{N, NC, T}) where {N, NC, T}
	    mask = bmask(hr.ibit)
	    step = 1<<(hr.ibit-1)
	    step_2 = step*2
	    nbit = nqubits(reg)
	    for j = 0:step_2:size(reg.state, 1)-step
	        for i = j+1:j+step
	            λ = bfloat(readbit(i-1, hr.cbits...), nbits=nbit-1)
	            if λ >= hr.C_value
	                u = rotmat(λ, hr.C_value)
	                YaoArrayRegister.u1rows!(state(reg), i, i+step, u...)
	            end
	        end
	    end
	    reg
	end

end

# ╔═╡ 077b5a14-2e88-4d24-bff3-b17335f00f20
function Solve(A,b,neig,scaling)

	# Prepare state
	U = exp(2*π*1im.*A);
	isValid(A,U,b);
	b_rangle = ArrayReg(b);
	eig_rangle = zero_state(neig);
	ancilla = zero_state(1);
	x = ArrayReg(join(b_rangle,eig_rangle,ancilla));

	# Build circuit
	s = size(U)[1];
	nqubits = log(2,s);
	ntotal = Int(nqubits+neig+1);
	sub_prob1 = PE(Int(neig),Int(nqubits),U);
	sub_prob2 = CRot{neig+1}([2:neig+1...,],1,scaling);
	HHL = chain(ntotal,subroutine(ntotal,sub_prob1,[2:ntotal...,]),subroutine(ntotal,sub_prob2,[1:(neig+1)...,]),subroutine(ntotal,sub_prob1',[2:ntotal...,]));


	# Apply
	x = x|>HHL;
	x = x|>focus!(1:(neig+1)...) |> select!(1) |> state |> vec;
	x = x./scaling;
	return x;
end

# ╔═╡ e1171d70-fa6b-42bf-ac75-44d1536e1a26
function PreProcessing(nqubits,dim)
	N = 2^nqubits;
	A = Matrix(Tridiagonal(-1*ones(N-1),2*ones(N),-1*ones(N-1)));
	E = Matrix(1.0I,N,N);
	if(dim==1)
		b = ones(ComplexF64,N);
	elseif(dim==2)
		A = kron(A,E) + kron(E,A);
		b = ones(ComplexF64,N*N);
	elseif(dim==3)
		A = kron(A,E,E) + kron(E,A,E) + kron(E,E,A);
		b = ones(ComplexF64,N*N*N);
	end
	b  = b ./norm(b);
	return A,b
end

# ╔═╡ 56a6dfc0-cc91-4605-b613-b2887c7feb6f
function PostProcessing(x_approx)
	real_x = real(x_approx);
	return real_x;
end

# ╔═╡ f0982801-2532-41e0-90be-bb3e8d9c6c58
function driver(nqubits,dim=1)
	nqubits = 5;
	neig = 10;
	A,b = PreProcessing(nqubits,dim);
	x_actual = A^-1*b;
	scaling = minimum(eigvals(A) .|> abs ) * 0.7;
	x = Solve(A,b,neig,scaling);
	x = PostProcessing(x);
	x_actual = PostProcessing(x_actual);
	return (x,x_actual)
end

# ╔═╡ Cell order:
# ╠═cf4de4d2-b32f-42b8-828e-aa9a15e0924c
# ╠═5fe8135e-7f83-4b20-b36e-b593153e333c
# ╠═3b41ca30-d458-11eb-1db5-151b58667f27
# ╠═5d8f93fa-2598-4ae1-8d67-9befe6385e7d
# ╠═caff44b6-72b8-4eca-8e79-6b47f163612e
# ╠═7de7d7af-66a6-477d-a43b-3b394bcdf3ac
# ╠═077b5a14-2e88-4d24-bff3-b17335f00f20
# ╠═e1171d70-fa6b-42bf-ac75-44d1536e1a26
# ╠═56a6dfc0-cc91-4605-b613-b2887c7feb6f
# ╠═f0982801-2532-41e0-90be-bb3e8d9c6c58
