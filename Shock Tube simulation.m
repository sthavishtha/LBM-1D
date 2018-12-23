%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MATLAB code to simulate 1D NSE in a shock tube
%   Assumed : 
%   Delta x (lattice distance) = Delta t (lattice time step) = 1 
%   c = 1, c_s (speed of sound) = 1/sqrt(3)
%   Periodic boundary conditions are applied at the corner grid points
%   Author : Sthavishtha Bhopalam Rajakumar
%   Updated date : 30-09-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear

%defined parameters in the question
nu = 0.06; %Diffusion coefficient
% u_max = 0; %macroscopic fluid velocity
n = 800; %number of grids
t_max = 500 ; %maximum time/iteration to be reached
c_s = 1./sqrt(3); %speed of sound
% ma= u_max./c_s; %mach number
tau = nu./(c_s*c_s); %relaxation time
beta = 1./(2*tau+1); %factor considering the effect of relaxation time
c = [0 1 -1]; %direction velocities

%initializing initial populations
rho = zeros(n,1); %particle densities
u = zeros(n,1); %particle velocities
for i=1:n %initial densities
    if i <= n/2
        rho(i) = 1.5;
    else
        rho(i) = 1.0;
    end
end
feq = zeros(3,n); %equilibrium populations
feq(1,:) = 2*rho(:)/3; %defining equilibrium populations
feq(2,:) = rho(:)./6;
feq(3,:) = rho(:)./6;
f = feq; %initial population is in equilbrium

for t = 1:t_max   %t = time/iteration
    
    if t~=1
        rho(:) = f(1,:) + f(2,:) + f(3,:); %density kinetic equation
        u(:) = (f(1,:)*c(1) + f(2,:)*c(2) + f(3,:)*c(3))./(rho(:))'; %momentum kinetic equation
        ma = u./(c_s); %local mach number
        feq(1,:) = 2*rho(:).*(1 - (u.^2/(2*c_s*c_s)))./3; %defining equilibrium populations - rewritten
        feq(2,:) = rho(:).*((1 + (u(:))./(c_s*c_s) + u.^2/(c_s*c_s)))./6;
        feq(3,:) = rho(:).*((1 - (u(:))./(c_s*c_s) + u.^2/(c_s*c_s)))./6;
%         feq(1,:) = 2*rho(:).*(2 - sqrt(1 + ma.^2))./3; %defining equilibrium populations
%         feq(2,:) = rho(:).*((-1./2 + u.*c(2)/(c_s*c_s) + sqrt(1 + ma.^2)))./3;
%         feq(3,:) = rho(:).*((-1./2 - u.*c(2)/(c_s*c_s) + sqrt(1 + ma.^2)))./3;
        f(1,:) = f(1,:) - 2*beta*(f(1,:) - feq(1,:)); %performing collision 
        f(2,:) = f(2,:) - 2*beta*(f(2,:) - feq(2,:));
        f(3,:) = f(3,:) - 2*beta*(f(3,:) - feq(3,:));
    end
    
    %advection across direction of c
    for i=n:-1:2
        f(2,i) = f(2,i-1);
    end
    
    for i=1:n-1
        f(3,i) = f(3,i+1);
    end
    
    %bounce-back
    f(2,1) = f(3,1);
    f(3,n) = f(2,n);
    
end
% for i=1:n
%     re = u.*i/nu;
% end
x = 1:n;
plot(x,rho(x),'k-');


