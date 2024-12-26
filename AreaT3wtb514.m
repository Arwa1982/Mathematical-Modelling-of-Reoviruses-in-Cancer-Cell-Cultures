load("T3wtb514.mat")
Tf = 120; % 120 hour
L = 30; % Upper boundary
M = 300; % M is the number of spaces between points a and b.
dx = 0.1; % (b-a)/M; % dx is delta x
dy = 0.1; % (b-a)/M; % dy is delta y
x = linspace(0, L, M+1); % M+1 equally spaced x vectors including a and b.
y = linspace(0, L, M+1);

% Time stepping

dt = 0.001; % (dx^2)/2D; % dt is delta t the time step
N = Tf/dt; % N is the number of time steps in the interval [0,30]
bn=0;

% Maximum ui and vn
maxun=max(un(:));
maxvn=max(vn(:));

% Compute the Area

 % Find the dynamic threshold based on the maximum value of ui
Threshold = 0.01;

% Find the radius and area of the circle with the dynamic threshold
for i = floor(M/2)+1:M+1
    for j = floor(M/2)+1:M+1
        if un(i, j) >= Threshold
            r1 = x(i);
            r2 = y(j);
            r3 = un(i, j);
            bn = 1;
            break
        end
    end
    if bn == 1
        break
    end
end
    

    % Compute the radius and area of the circle
    if bn == 1
        radius = sqrt((r1 - 15)^2 + (r2 - 15)^2);
        areacircle = pi * radius^2;
    end

 % Plot results
figure;
pn=pcolor(x, y,un);
set(pn, 'EdgeColor', 'none');
xlabel('x');
ylabel('y');
colorbar;
grid on;
set(gca,'FontSize',24);
hold on;  % Enable the "hold" to overlay additional graphics on the same figure
% Plot the circle
rectangle('Position', [15-radius, 15-radius, 2*radius, 2*radius], ...
          'Curvature', [1, 1], 'EdgeColor', 'r','LineWidth', 2.5);
% Plot a dummy plot for the legend
plot(NaN, NaN, 'r', 'LineWidth', 2, 'DisplayName', 'Selected Area');
hold off;  % Release the "hold" to avoid overlaying on future plots


figure;
pn=pcolor(x, y,vn);
set(pn, 'EdgeColor', 'none');
xlabel('x');
ylabel('y');
colorbar;
grid on;
set(gca,'FontSize',24);
hold on;  % Enable the "hold" to overlay additional graphics on the same figure
% Plot the circle
rectangle('Position', [15-radius, 15-radius, 2*radius, 2*radius], ...
          'Curvature', [1, 1], 'EdgeColor', 'r','LineWidth', 2.5);
% Plot a dummy plot for the legend
plot(NaN, NaN, 'r', 'LineWidth', 2, 'DisplayName', 'Selected Area');
hold off;  % Release the "hold" to avoid overlaying on future plots


% Initialize variables to store coordinates
x_coordinates = [];
y_coordinates = [];

% Find all grid points with coordinates (r1, r2)
for i = floor(M/2)+1:M+1
    for j = floor(M/2)+1:M+1
        if abs(x(i) - r1) < 0.01 && abs(y(j) - r2) < 0.01
            x_coordinates = [x_coordinates, x(i)];
            y_coordinates = [y_coordinates, y(j)];
        end
    end
end

% Compute vn at the found coordinates and sum them up
total_vn = 0;
for k = 1:length(x_coordinates)
    % Find the indices corresponding to the coordinates
    [~, i_index] = min(abs(x - x_coordinates(k)));
    [~, j_index] = min(abs(y - y_coordinates(k)));
    
    % Sum the vn values at those indices
    total_vn = total_vn + vn(i_index, j_index);
end

% Display the total vn value
disp('Total vn value at radius 2.2:');
disp(total_vn);


figure;
s=surf(x,y,vn)
s.EdgeColor = 'none';
colorbar;
grid on;
set(gca,'FontSize',24);