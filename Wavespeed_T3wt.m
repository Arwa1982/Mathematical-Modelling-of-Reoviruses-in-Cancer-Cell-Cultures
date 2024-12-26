% Solve the PDE system of virus-infected cells system for T3wt
% Gamma=0.96 and b=514
% Analytic wave speed= 0.0437 and numerical wave speed=0.0437
% Space and time
x = linspace(0,20,2001); % dx=0.01,  
t = linspace(0,458,91601); % dt=0.005

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
plot(x,V(i,:),'Color','[0.6350 0.0780 0.1840]','LineWidth',2)
end
hold off
ylim([0 1])
xlabel('x (mm)')
ylabel('V_{T3wt}(x,t)')
set(gca,'FontSize',24)
grid on



VT = V(1: 14400:length(t),:);
g =[];
for i = 1:7
linearIndices = find(VT(i,:)<0.001);
g(i) = linearIndices(1);
end

xx = [];
for i = 1:7
xx(i) = x(g(i));
end
TT = t(1: 14400:length(t));    
sp = [];
for i = 1:7
    if i<7
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
b=514;
gamma=0.96;
nu=0.01;
kc=1/4.14;
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
  u0(2)=0.1684;
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