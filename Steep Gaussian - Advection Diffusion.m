%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MATLAB code to simulate 1D advection-diffusion for a steep gaussian profile
%   Assumed : 
%   Delta x (lattice distance) = Delta t (lattice time step) = 1 
%   c = 1, c_s (speed of sound) = 1/sqrt(3)
%   Periodic boundary conditions are applied at the corner grid points
%   Author : Sthavishtha Bhopalam Rajakumar
%   Updated date : 30-09-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear

%defined parameters in the question
d = 5e-8; %Diffusion coefficient
u_max = 0.1; %macroscopic fluid velocity
n = 1000; %number of grids
t_max = 4000; %maximum time/iteration to be reached
c_s = 1./sqrt(3); %speed of sound
ma = u_max./c_s; %mach number
tau = d./(c_s*c_s); %relaxation time
beta = 1./(2*tau+1); %factor considering the effect of relaxation time

%initializing initial populations
rho = zeros(n,1); %particle densities
for i=1:n
    rho(i) = 1.0+0.5*exp(-5000*power(((i-1)/(n-1) - 0.25),2)); %initial density variation
end
feq = zeros(3,n); %equilibrium populations
feq(1,:) = 2*rho(:).*(2 - sqrt(1 + ma*ma))./3; %defining equilibrium populations
feq(2,:) = rho(:).*((u_max - c_s*c_s)/(2*c_s*c_s) + sqrt(1 + ma*ma))./3;
feq(3,:) = rho(:).*((-u_max - c_s*c_s)/(2*c_s*c_s) + sqrt(1 + ma*ma))./3;
f = feq; %initial population is in equilibrium 

for t = 1:t_max     %t = time/iteration
    
    if t~=1
        rho(:) = f(1,:) + f(2,:) + f(3,:); %density kinetic equation
        feq(1,:) = 2*rho(:).*(2 - sqrt(1 + ma*ma))./3; %defining equilibrium populations
        feq(2,:) = rho(:).*((u_max - c_s*c_s)/(2*c_s*c_s) + sqrt(1 + ma*ma))./3;
        feq(3,:) = rho(:).*((-u_max - c_s*c_s)/(2*c_s*c_s) + sqrt(1 + ma*ma))./3;
        f(1,:) = f(1,:) - 2*beta*(f(1,:) - feq(1,:)); %performing collision 
        f(2,:) = f(2,:) - 2*beta*(f(2,:) - feq(2,:));
        f(3,:) = f(3,:) - 2*beta*(f(3,:) - feq(3,:));
    end
   
    temp_left = f(3,1); %periodic boundary conditions
    temp_right = f(2,n);
    
    %advection across direction of c
    for i=n:-1:2
        f(2,i) = f(2,i-1);
    end
    
    for i=1:n-1
        f(3,i) = f(3,i+1);
    end
    
    f(2,1) = temp_right;
    f(3,n) = temp_left;
    
end
%postprocessing
x = 1:n;
rho_act(x) = 1+0.5*exp(-5000*power(((x-1)/(n-1) - 0.25),2));
plot(x,rho_act(x),'k--',x,rho(x),'k-');

