clear all
close all

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

% Constant Values
D = 0.01;
gamma = 0.96;
nu = 0.01;
alpha = 0.057;
b=514;

   % Pre-allocation
    unp1 = zeros(M+1, M+1); % C(x,t)
    uip1 = zeros(M+1, M+1); % I(x,t)
    vnp1 = zeros(M+1, M+1); % V(x,t)

   
   % Initial Conditions
   un(1:M+1, 1:M+1) = 1; % IC's of C(x,t)
   ui(1:M+1, 1:M+1) = 0; % IC's of I(x,t)
   mu = [15  15];
   Sigma = [0.01 0; 0 0.01];
  [X1, X2] = meshgrid(x, y);
   X = [X1(:) X2(:)];
  vn = mvnpdf(X, mu, Sigma);
  vn = reshape(vn, length(y), length(x));



   for n = 1:N
    % Boundary conditions on left, right, bottom, and top flux are zero
     un(1, 1:M+1) = un(2, 1:M+1); % Boundary conditions on left flux is zero  
     ui(1, 1:M+1) = ui(2, 1:M+1); % Boundary conditions on left flux is zero  
     vn(1, 1:M+1) = vn(2, 1:M+1); % Boundary conditions on left flux is zero  

     un(M+1, 1:M+1) = un(M, 1:M+1); % Boundary conditions on right flux is zero  
     ui(M+1, 1:M+1) = ui(M, 1:M+1); % Boundary conditions on right flux is zero  
     vn(M+1, 1:M+1) = vn(M, 1:M+1); % Boundary conditions on right flux is zero  

     un(1:M+1, 1) = un(1:M+1, 2); % Boundary conditions on bottom flux is zero  
     ui(1:M+1, 1) = ui(1:M+1, 2); % Boundary conditions on bottom flux is zero  
     vn(1:M+1, 1) = vn(1:M+1, 2); % Boundary conditions on bottom flux is zero  

     un(1:M+1, M+1) = un(1:M+1, M); % Boundary conditions on top flux is zero  
     ui(1:M+1, M+1) = ui(1:M+1, M); % Boundary conditions on top flux is zero  
     vn(1:M+1, M+1) = vn(1:M+1, M); % Boundary conditions on top flux is zero  


    for i = 2:M
        for j = 2:M

            % Source function for u and v
            srcu = -gamma*nu*vn(i, j)*un(i, j);
            srci = gamma*nu*vn(i, j)*un(i, j)-alpha*ui(i, j);
            srcv = alpha*b*ui(i, j)-gamma*vn(i, j);

            
            vxx = (vn(i-1, j)-2*vn(i, j)+vn(i+1, j))/dx^2; % Laplacian v
           
            vyy = (vn(i, j-1)-2*vn(i, j)+vn(i, j+1))/dy^2; % Laplacian v
            
            Lapv = vxx + vyy;

            unp1(i, j) = un(i, j) + dt*srcu; % Full dynamics 
            uip1(i, j) = ui(i, j) + dt* srci; % Full dynamics
            vnp1(i, j) = vn(i, j) + dt*(D*Lapv + srcv); % Full dynamics
        end
    end

    un = unp1; % Updated values
    ui = uip1; % Updated values
    vn = vnp1; % Updated values
end
           

    % Compute plaque size for the current b value
    TotalC = sum(un(:));
    TotalI = sum(ui(:));
    TotalV = sum(vn(:));
    maxui= max(ui(:));
     maxvn= max(vn(:));
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


figure;
pui=pcolor(x, y,ui);
set(pui, 'EdgeColor', 'none');
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
pun=pcolor(x, y,vn);
set(pun, 'EdgeColor', 'none');
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

