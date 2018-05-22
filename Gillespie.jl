Pkg.add("Plots")
using Plots;

plotly();


for Thorizon=1:20
	rng = MersenneTwister(123445454444888985559541215652);
	X = [0.0, 0.0]


	E=5;
	S=500;

	Kd=10;
	kdg=2;

	r=0.4;
	kappa1=5;
	kappa2=1;
	b=kdg/(E*S);

	v=1;
	u=v*b*S*S;

	k1=kappa1*S*S*b;
	k2=kappa2*b*E*S;


	R=0.4/(kdg*sqrt(Kd));

	Rh=R*b*E*S*S;


	function Rates(X, Rh, k1, k2, b, u,E)
		W=[0.0;0.0;0.0;0.0];
		W[1]=Rh+k1*X[2];
		W[2]=k2*X[1];
		W[3]=b*X[1]*(X[1]-1)*(E-X[2]);
		W[4]=u*X[2];
		return W;
	end


	W=[0.0; 0.0; 0.0 ;0.0];
	C=[1.0 0.0 0.0 0.0 ; 1.0 1.0 0.0 0.0; 1.0 1.0 1.0 0.0; 1.0 1.0 1.0 1.0];
	rvec=Vector{Float64}[[1,0], [-1, 0], [-2, 1], [2, -1]];
	timevec=Float64[];
	Xvec=Vector{Float64}[X];
	tau=5;
	t=0; 	
	n=1;
	push!(timevec, t);

	while t<Thorizon
		W=Rates(Xvec[n], Rh, k1, k2, b, u, E);
		Wj=*(C, W);
		W0=Wj[4];

		X = rand(rng, 2);

		tau=(1/W0)*log(1/(X[1]));

		j=1;

		while (j<length(Wj)) && (X[2]*W0>=Wj[j])
			j+=1;
		end 
		t=t+tau;
		Xtmp=Xvec[n]+rvec[j];
		push!(timevec, t);
		push!(Xvec, Xtmp);
		n=n+1;
	end
	print(Thorizon);
	print(",");	
	println(n);
end
X1G=Float64[]
X2G=Float64[]
for i=1:length(Xvec)
	push!(X1G, Xvec[i][1]);
	push!(X2G, Xvec[i][2]);
end


plot(timevec, X1G, label="Fast Rate Variable")
plot(timevec, X2G, label="Slow Rate Variable")

plot(timevec, X1G, label="Fast Rate Variable")
plot!(timevec, X2G, label="Slow Rate Variable")

boxplot(X1G, lab="X1")
boxplot!(X2G, lab="X2")

boxplot(X1G, lab="Gillespie")
boxplot!(X1, lab="QSS Gillespie")

boxplot(X2G, lab="Gillespie")
boxplot!(X2, lab="QSS Gillespie")
