function [cp, ctau, cd, cl] = coeff_maxwell(param_eq,  delta)

alpha = param_eq.alpha;
Tw = param_eq.Tw;
Tinf = param_eq.Tinf;
s = param_eq.s;

f = 1-alpha;
theta = pi/2-delta;

cd = 2*((1-f*cos(2*theta))./(sqrt(pi)*s)).*exp(-s^2*sin(theta).^2)...
    + (sin(theta)/s^2).*(1+2*s^2+f*(1-2*s^2*cos(2*theta))).*erf(s*sin(theta))...
    + (1-f)/s * sqrt(pi)*sin(theta).^2*sqrt(Tw/Tinf);

cl = ((4*f)/(sqrt(pi)*s))*sin(theta).*cos(theta).*exp(-s^2*sin(theta).^2)...
    + (cos(theta)/s^2).*(1+f*(1+4*s^2*sin(theta).^2)).*erf(s*sin(theta))...
    + ((1-f)/s)*sqrt(pi)*sin(theta).*cos(theta)*sqrt(Tw/Tinf);

cp   = cd.*cos(delta) - cl.*sin(delta);
ctau = cd.*sin(delta) + cl.*cos(delta);