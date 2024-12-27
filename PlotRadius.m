% Parameter values
D=0.01;
t=1;
chi=1000;
gamma=0:0.1:15;


% Function
 f=@(gamma) sqrt(4*D*t*(-log(4*pi*D*t/chi)-gamma*t));


 % Simulation
 figure()
 fplot(f,'k','LineWidth',2.5)
 hold on
 xline(0.28,'r','LineWidth',2.5)
 xline(0.96,'b','LineWidth',2.5)
 hold off
 xlim([0,10])
 ylim([0,0.7])
 xlabel('$\gamma_{b}$','Interpreter','latex')
 ylabel('r (mm)')
 set(gca,'FontSize',24)
 grid on




