% Parameters values
d = 0.01;
a = 0.057;
v = 0.01;

% Define the range for sensitivity analysis
b_range = 500:25:1000;  % Range of values for parameter b
gamma_range = 0.13:0.001:0.5; % Range of values for parameter gamma
min_c_values = zeros(length(b_range), length(gamma_range)); % Initialize array to store minimum c values

% Perform sensitivity analysis
for i = 1:length(b_range)
    for j = 1:length(gamma_range)
        b = b_range(i);
        gamma = gamma_range(j);
        
        % Define the function for SV5 virus
        fun = @(x) (d * x.^2 * (a + x)) ./ (x.^2 + (a + gamma) * x - a * gamma * (b * v - 1));
        % Set options for optimization
        opts = optimset('MaxFunEvals', 1e15); % Increase the maximum number of function evaluations
        % Find minimum
        vals = fminsearch(fun, 0.3); % Ensure positive and real
        minresult = fun(vals);
        if isreal(minresult) && ~isnan(minresult) && minresult > 0 % Check if the minimum value is real, positive, and not NaN
            min_c_values(i, j) = sqrt(minresult);
        else
            min_c_values(i, j) = NaN; % If not real, positive, or NaN, mark it as NaN
        end
    end
end

% Plot results

figure(1)
surf(gamma_range, b_range, min_c_values,'edgecolor', 'none')
xlabel('\gamma_b')
ylabel('$\tilde{b}$','Interpreter','latex')
zlabel('c^{*}')
set(gca,'Fontsize',24)



figure(2)
surf(gamma_range, b_range, min_c_values,'edgecolor', 'none')
xlabel('\gamma_b')
ylabel('$\tilde{b}$','Interpreter','latex')
zlabel('c^{*}')
set(gca,'Fontsize',24)
grid on

