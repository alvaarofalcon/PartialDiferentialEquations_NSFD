clear all;
close all;
clc;

tic
tiempo_a_graficar = 1000;
save_interval = 1000;

params.lambda_S = 0.1;params.lambda_E = 0.1; params.lambda_I = 0.1; 
params.lambda_Q = 0.1;params.lambda_R = 0.1; params.beta1  = 0.001;
params.beta2 = 0.003; params.mu = 0.06; params.delta = 0.01; 
params.gamma = 0.04; params.alfa = 0.04; params.rho = 0.02;
params.Lambda = 1; params.h = 0.01; params.k = 0.0001; %pasos espacio y tiempo 
% respectivamente

total_pasos = tiempo_a_graficar / params.k;
t_exp = 1:total_pasos;

%% Espacio 1D discretizado

S0 = 0.15 * ones(1, 101);
E0 = 0.01 * ones(1, 101);
I0 = 0.005 * ones(1, 101);
Q0 = zeros(1, 101);
R0 = zeros(1, 101);

%% Vector tiempo
[S_hist, E_hist, I_hist, Q_hist, R_hist] = OneDScheme(t_exp, params, S0, E0, I0, Q0, R0,...
    save_interval);
%% Grafico de malla en 1D

Nx = length(S0);                 % número de puntos espaciales (101)
Nt = length(S_hist);             % número de pasos guardados
x = linspace(0, 1, Nx);          % dominio espacial, ajústalo a lo que uses
t = linspace(0, tiempo_a_graficar, Nt);  % eje temporal

% Convertir cell array a matriz 2D: tiempo × espacio
S_mat = cell2mat(S_hist');  
E_mat = cell2mat(E_hist');
I_mat = cell2mat(I_hist');
Q_mat = cell2mat(Q_hist');
R_mat = cell2mat(R_hist');

[X, T] = meshgrid(x, t);

%% Graficar
figure('Name', 'Evolución 1D en el tiempo');
tiledlayout(3, 2, 'TileSpacing', 'compact');
colormap("jet");

nexttile;
surf(X, T, S_mat, 'EdgeColor', 'none');
set(gca, 'XDir', 'reverse');
xlabel('Space'); ylabel('Time'); zlabel('S');
title('Susceptible'); view(45,60); shading interp;
clim([0 1]);

nexttile;
surf(X, T, E_mat, 'EdgeColor', 'none');
set(gca, 'XDir', 'reverse');
xlabel('Space'); ylabel('Time'); zlabel('E');
title('Exposed'); view(45,60); shading interp;
clim([0 1]);

nexttile;
surf(X, T, I_mat, 'EdgeColor', 'none');
set(gca, 'XDir', 'reverse');
xlabel('Space'); ylabel('Time'); zlabel('I');
title('Infected'); view(45,60); shading interp;
clim([0 1]);

nexttile;
surf(X, T, Q_mat, 'EdgeColor', 'none');
set(gca, 'XDir', 'reverse');
xlabel('Space'); ylabel('Time'); zlabel('Q');
title('Quarantined'); view(45,60); shading interp;
clim([0 1]);

nexttile;
surf(X, T, R_mat, 'EdgeColor', 'none');
set(gca, 'XDir', 'reverse');
xlabel('Space'); ylabel('Time'); zlabel('R');
title('Recovered'); view(45,60); shading interp;
clim([0 1]);

cb = colorbar;
cb.Layout.Tile = 'east';