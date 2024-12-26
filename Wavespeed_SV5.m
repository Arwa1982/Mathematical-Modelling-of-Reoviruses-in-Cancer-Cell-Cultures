% Solve the PDE system of virus-infected cells system for SV5
% Gamma=0.28 and b=732
% Analytic wave speed= 0.0594 and numerical wave speed=0.0599
% Space and time
x = linspace(0,20,2001); % dx=0.01,  
t = linspace(0,337,67401); % dt=0.005

% PDEPE Solver
m = 0;
sol = pdepe(m,@angiopde,@angioic,@angiobc,x,t);
% Solution
V = sol(:,:,1);
I = sol(:,:,2);

% Plot and simulation

figure(1)
hold on
for i=1: 14400:length(t) % 3 days
plot(x,V(i,:),'Color','[0.4940 0.1840 0.5560]','LineWidth',2)
end
hold off
ylim([0 1])
xlabel('x (mm)')
ylabel('V_{SV5}(x,t)')
set(gca,'FontSize',24)
grid on


VT = V(1: 14400:length(t),:);
g =[];
for i = 1:5
linearIndices = find(VT(i,:)<0.001);
g(i) = linearIndices(1);
end

xx = [];
for i = 1:5
xx(i) = x(g(i));
end
TT = t(1: 14400:length(t));    
sp = [];
for i = 1:5
    if i<5
    sp(i) = ((xx(i+1)-xx(i))/(TT(i+1)-TT(i)));
    end
end

format long
mean(sp)
% ------------------------------------------------
function [c,f,s] = angiopde(x,t,u,dudx) % Equation to solve
dv=0.01;
di=0;
a=0.057;
b=732;
gamma=0.28;
nu=0.01;
kc=1/6.32;
c = [1; 1];
f = [dv; di] .* dudx;
s = [a*b*u(2)-gamma*u(1)-gamma*u(1).^2./kc;
     gamma*nu*u(1)-a*u(2)];
end
% ---------------------------------------------
function u0 = angioic(x) % Initial Conditions
u0 = [0; 0];
if x >= 0 && x <=0.5
  u0(1)=1;
  u0(2)=0.04912;
end
end
% ---------------------------------------------
function [pl,ql,pr,qr] = angiobc(xl,ul,xr,ur,t) % Boundary Conditions
pl = [0; 0];
ql = [1; 1];
pr = pl;
qr = ql;
end
% ---------------------------------------------