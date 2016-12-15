function modalparam_test()
	# State space model for testing modalparam 
	# and modalparami
	# 
	# javier.cara@upm.es

	M = [35.0 0.0;0.0 17.5] # mass matrix
	K = [12250.0 -3500.0;-3500.0 3500.0] # stiffness matrix
	g = [0.02;0.02] # damping ratios
	
	# teorethical eigenvalues and eigenvectors
	W2,V = eig(K,M)
	W = sqrt(W2)
	
	# damping matrix
	Mm = V'*M*V
	Gm = 2*Mm*diagm(W)*diagm(g)
	G = inv(V')*Gm*inv(V) # damping matrix
	
	# state-space matrices
	Minv = [ 1/M[1,1] 0.0;0.0 1/M[2,2] ]
	Ac = [zeros(2,2) eye(2,2);-Minv*K -Minv*G] # continuous A matrix
	dt = 0.02 # time step
	A = expm(Ac*dt) # discrete A matrix
	C = [-Minv*K -Minv*G] # discrete C matrix
	Bc = [zeros(2,2);Minv]
	B = (A - eye(4))*inv(Ac)*Bc
	
	# eigenvectors with max. component = 1
	V1 = zeros(2,2)
	for j in 1:2
		maxval = V[1,j]
		if abs(V[2,j]) > abs(maxval)
			maxval = V[2,j]
		end
		V1[:,j] = V[:,j]/maxval
	end	
	# updated mass masses
	Mm1 = V1'*M*V1
	mm1 = diag(Mm1)
	
	ssm = Dict("M"=>M,"K"=>K,"G"=>G,"wm"=>W,"zm"=>g,"Vm"=>V1,"mm"=>mm1,
							"dt"=>dt,"A"=>A,"B"=>B,"C"=>C)
	
	return ssm
	
end

