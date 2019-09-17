function [] = test_simulation()

	%% Cauchy parameters
	a = [0.5 0];
	b = [1 1];
	p = 0.8;

	%% Triangle vertices
	x1 = [-2 -2];
	x2 = [0 4];
	x3 = [4 0];
	
	%% plot axes
	AXES = [-5 5 -5 5]	
	PLOT_RESOLUTION = 64;
	
	%% plot the Cauchy distribution inside the triangle
	distrib_reference = zeros(PLOT_RESOLUTION, PLOT_RESOLUTION);
	for j=1:PLOT_RESOLUTION
	for i=1:PLOT_RESOLUTION

		x(1) = AXES(1) + (AXES(2)-AXES(1)) * (i-0.5)/PLOT_RESOLUTION;
		x(2) = AXES(1) + (AXES(2)-AXES(1)) * (j-0.5)/PLOT_RESOLUTION;
		distrib_reference(i,j) = Cauchy(x, a, b, p) * InsideTriangle(x1, x2, x3, x);
			
	end
	end
	figure
	imagesc(distrib_reference)
	
	%% generate 8k random points and plot
	numberOfPoints = 8192
	distrib_simulated = zeros(PLOT_RESOLUTION, PLOT_RESOLUTION);
	for i = 1:numberOfPoints
			
		% generate random numbers
		U1 = rand;
		U2 = rand;
		
		% get point
		x = SimulateCauchyOverPolygon_analytic(x1, x2, x3, a, b, p, U1, U2);

		% splat point in distribution
		i = floor((x(1)-AXES(1)) / (AXES(2)-AXES(1)) * PLOT_RESOLUTION) + 1;
		j = floor((x(2)-AXES(3)) / (AXES(4)-AXES(3)) * PLOT_RESOLUTION) + 1;
		if i>=1 && i <= PLOT_RESOLUTION && j>=1 && j <= PLOT_RESOLUTION
			distrib_simulated(i,j) += 1;
		endif
	
	end
	figure
	imagesc(distrib_simulated)
	
endfunction

function [inside] = InsideTriangle(x1, x2, x3, x)

	x12 = x2-x1;
	x23 = x3-x2;
	x31 = x1-x3;
	x10 = x-x1;
	x20 = x-x2;
	x30 = x-x3;
	
	inside = true;
	inside = inside & (-sign(x12(1)*x31(2)-x12(2)*x31(1)) == sign(x12(1)*x10(2)-x12(2)*x10(1)));
	inside = inside & (-sign(x23(1)*x12(2)-x23(2)*x12(1)) == sign(x23(1)*x20(2)-x23(2)*x20(1)));
	inside = inside & (-sign(x31(1)*x23(2)-x31(2)*x23(1)) == sign(x31(1)*x30(2)-x31(2)*x30(1)));

endfunction


function [x] = SimulateCauchyOverPolygon_analytic(x1, x2, x3, a, b, p, U1, U2)

	%% 
	x1_std = CauchyToStd(x1, a, b, p);
	x2_std = CauchyToStd(x2, a, b, p);
	x3_std = CauchyToStd(x3, a, b, p);
	
	%% Spherical triangle vertices
	w1 = x1_std / norm(x1_std);
	w2 = x2_std / norm(x2_std);
	w3 = x3_std / norm(x3_std);
	
	%% Simulate point uniformly in the spherical triangle 
	w = SimulateSphericalTriangle(w1, w2, w3, U1, U2);
	
	%% Map point from hemisphere to plane
	x_std = w / w(3);
	
	%% 
	x = StdToCauchy(x_std, a, b, p);

endfunction

%% Evaluate the Cauchy distribution
function [v] = Cauchy(x, a, b, p)
	
	x_std = CauchyToStd(x, a, b, p);
	v = CauchyStd(x_std) * 1 / (b(1)*b(2)*sqrt(1-p^2));
	
endfunction

%% Evaluate the standard Cauchy distribution
function [v] = CauchyStd(x)

	v = 1/(2*pi) * 1/(x(1)^2+x(2)^2+1)^(3/2);
	
endfunction

%% Location-scale correlation forward
function [x_std] = CauchyToStd(x_c, a, b, p)

	x_std(1) = (x_c(1)-a(1)) / b(1);
	x_std(2) = (b(1)*(x_c(2)-a(2)-p*b(2)*(x_c(1)-a(1)))) / (b(1)*b(2)*sqrt(1-p^2));
	x_std(3) = 1;
	
endfunction

%% Location-scale correlation backward
function [x_c] = StdToCauchy(x_std, a, b, p)

	x_c(1) = b(1)*x_std(1) + a(1);
	x_c(2) = b(2) * (p*x_std(1) + x_std(2)*sqrt(1-p^2)) + a(2);
	x_c(3) = 1;

endfunction

%% Generate a random point uniformly in a spherical triangle 
%% Based on "Stratified Sampling of Spherical Triangles" by James Arvo
function [w] = SimulateSphericalTriangle(w1, w2, w3, U1, U2)

	%%
	a = acos(dot(w2, w3));
	b = acos(dot(w1, w3));
	c = acos(dot(w1, w2));

	%% Spherical angles of the triangle
	alpha = SphericalAngle(w1, w2, w3);
	beta = SphericalAngle(w2, w1, w3);
	gamma = SphericalAngle(w3, w1, w2);

	%% solid angle 
	Omega = U1 * (alpha + beta + gamma - pi);
	
	%% 
	s = sin(Omega - alpha);
	t = cos(Omega - alpha);
	u = t - cos(alpha);
	v = s + sin(alpha)*cos(c);
	q = ((v*t-u*s) * cos(alpha) - v) / ((v*s+u*t)*sin(alpha));
	C2 = q*w1 + sqrt(1-q*q) * normalize(w3 - dot(w3,w1)*w1);
	z = 1-U2*(1-dot(C2, w2));
	
	%% result
	w = z*w2 + sqrt(1-z*z) * normalize(C2 - dot(C2,w2)*w2);
	
endfunction

%% Angle at w1 of the spherical triangle w1, w2, w3 
function [alpha] = SphericalAngle(w1, w2, w3)

	w12 = normalize((w2-w1) - w1*dot(w2-w1, w1));
	w13 = normalize((w3-w1) - w1*dot(w3-w1, w1));
	alpha = acos(dot(w12, w13));
	
endfunction

%% Normalize 3d vector
function [w] = normalize(w)

	w = w / norm(w);

endfunction
