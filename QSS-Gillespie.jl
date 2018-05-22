Pkg.add("Plots")
using Plots;

plotly();

for Thorizon=1:20
	rng = MersenneTwister(123445454444888985559541215652);
	X = 0.0

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

	function Rates(X, Rh, k1, k2, b, u, E, S, v)
		W=[0.0; 0.0; 0.0; 0.0];

		X2=E*((X/S)^2);
		X2=X2/(v+(X/S)^2);

		W[1]=Rh+k1*X2;
		W[2]=k2*X;
		W[3]=b*X*(X-1)*(E-X2);
		W[4]=u*X2;

		return W;
	end


	W=[0.0; 0.0; 0.0 ;0.0];
	C=[1.0 0.0 0.0 0.0 ; 1.0 1.0 0.0 0.0; 1.0 1.0 1.0 0.0; 1.0 1.0 1.0 1.0];

	rvec=Float64[1; -1; -2; 2];
	timevec=Float64[];
	Xvec=Float64[X];
	tau=0;
	t=0; 	
	n=1;
	push!(timevec, t);

	while t<Thorizon
		W=Rates(Xvec[n], Rh, k1, k2, b, u, E, S, v);
		Wj=*(C, W);
		W0=Wj[4];

		Xran = rand(rng, 2);

		tau=(1/W0)*log(1/(Xran[1]));

		j=1;

		while (j<length(Wj)) && (Xran[2]*W0>=Wj[j])
			j+=1;
		end 

		t=t+tau;

		Xtmp=Xvec[n]+rvec[j];
		push!(timevec, t);
		push!(Xvec, Xtmp);
		n=n+1;
	#Change the value of X2 at n iteration here
	end

	print(Thorizon);
	print(",");	
	println(n);
end

X1=Float64[];
X2=Float64[];
for i=1:length(Xvec)
	push!(X1, Xvec[i]);
	push!(X2, E*Xvec[i]*Xvec[i]/(v+Xvec[i]*Xvec[i]));
end


plot(timevec, X1, label="Slow Rate Variable")
plot!(timevec, X2, label="Fast Rate Variable")
plot(timevec, X1, label="Slow Rate Variable")
plot(timevec, X2, label="Fast Rate Variable", ylimits=(0,6))
boxplot(X1, lab="X1")
boxplot!(X2, lab="X2")
