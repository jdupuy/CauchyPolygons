%% Compares the analytic polygonal integral of a bivariate Cauchy distribution
%% against a Monte Carlo integration
function [] = test_integration()

	%% Cauchy parameters
	a = [2 0.3];
	b = [0.3 2];
	p = 0.5;

	%% Triangle vertices
	x1 = [0 0];
	x2 = [1 0];
	x3 = [0 1];
 
	%% Analytic integration
	I_analytic = IntegrateCauchyOverPolygon_analytic(x1, x2, x3, a, b, p)

	%% Numerical integration (Monte Carlo)	
	NumberOfSamples = 10000;
	I_numeric = IntegrateCauchyOverPolygon_numeric(x1, x2, x3, a, b, p, NumberOfSamples)

endfunction

function [I] = IntegrateCauchyOverPolygon_analytic(x1, x2, x3, a, b, p)

	%% 
	x1_std = CauchyToStd(x1, a, b, p);
	x2_std = CauchyToStd(x2, a, b, p);
	x3_std = CauchyToStd(x3, a, b, p);
	
	%% spherical triangle vertices
	w1 = x1_std / norm(x1_std);
	w2 = x2_std / norm(x2_std);
	w3 = x3_std / norm(x3_std);
	
	%% Solid angle
	Omega = SolidAngleSphericalTriangle(w1, w2, w3);
	
	%% Normalize
	I = Omega / (2*pi);	

endfunction

function [I] = IntegrateCauchyOverPolygon_numeric(x1, x2, x3, a, b, p, NumberOfSamples)

	%% triangle area
	P1 = [x1(1) x1(2) 1];
	P2 = [x2(1) x2(2) 1];
	P3 = [x3(1) x3(2) 1];
	Area = 0.5 * norm(cross(P2-P1, P3-P1));

	%% MC integration
	sum = 0;
	for sample = 1:NumberOfSamples

		%% generate point x uniformly inside the triangle
		u1 = rand;
		u2 = rand;
		x = (1-sqrt(u1))*x1 + sqrt(u1)*(u2*x2 + (1-u2)*x3);
		
		%% eval Cauchy at x
		sum += Cauchy(x, a, b, p);
			
	end

	%% display	
	I = Area * sum/NumberOfSamples;
		
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

%% Normalize 3d vector
function [w] = normalize(w)

	w = w / norm(w);

endfunction

%% Solid angle of a spherical triangle 
function [Omega] = SolidAngleSphericalTriangle(w1, w2, w3)

	%% Spherical angles of the triangle
	alpha1 = SphericalAngle(w1, w2, w3);
	alpha2 = SphericalAngle(w2, w1, w3);
	alpha3 = SphericalAngle(w3, w1, w2);
	
	%% solid angle 
	Omega = alpha1 + alpha2 + alpha3 - pi;

endfunction

%% Angle at w1 of the spherical triangle w1, w2, w3 
function [alpha] = SphericalAngle(w1, w2, w3)

	w12 = normalize((w2-w1) - w1*dot(w2-w1, w1));
	w13 = normalize((w3-w1) - w1*dot(w3-w1, w1));
	alpha = acos(dot(w12, w13));
	
endfunction