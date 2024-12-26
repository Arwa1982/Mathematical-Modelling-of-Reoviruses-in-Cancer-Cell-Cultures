% Finding the minimum wave speed with T3wt virus

% Parameters values
d = 0.01;
a = 0.057;
gammat = 0.96;
bt = 514;
v = 0.01;
gammas = 0.28;
bs = 732;

% Define the function for T3wt virus

funt = @(x) d * x.^2 * (a + x) / (x.^2 + (a + gammat) * x - a * gammat * (bt * v - 1));
valt=fminsearch(funt,0.3)
minT3wt=funt(valt)
cT3wt=sqrt(minT3wt)

% Define the function for SV5 virus
funs = @(x) d * x.^2 * (a + x) / (x.^2 + (a + gammas) * x - a * gammas * (bs * v - 1));
vals=fminsearch(funs,0.3)
minSV5=funs(vals)
cSV5=sqrt(minSV5)

% Simulation
figure

% Plot the T3wt virus function (ft)
fplot(funt, [0 0.5], 'r', 'LineWidth', 2.5)
hold on

% Plot the SV5 virus function (fs)
fplot(funs, [0 0.5],'Color','[0.4940 0.1840 0.5560]', 'LineWidth', 2.5)

hold off
xlim([0 0.5]);
ylim([-0.1 0.1]);
xlabel('$\varrho$','Interpreter','latex');
ylabel('$\chi(\varrho)$','Interpreter','latex')
grid on
%legend('T3wt Virus', 'SV5 Virus')  % Add legend to differentiate between the two functions
set(gca, 'fontsize', 24)

% Mark the minimum point for SV5 virus on the plot with a filled marker
hold on
plot(valt, minT3wt, 'r.', 'MarkerSize', 50, 'MarkerFaceColor', 'r', 'LineWidth', 2.5);
plot(vals, minSV5,'.','Color','[0.4940 0.1840 0.5560]', 'MarkerSize', 50, 'MarkerFaceColor','[0.4940 0.1840 0.5560]','LineWidth', 2.5);
hold off
legend('T3wt','SV5','$\varrho_{*T3wt}$','$\varrho_{*SV5}$','Interpreter','latex')