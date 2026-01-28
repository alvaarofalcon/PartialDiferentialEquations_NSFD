%% Parámetros iniciales
clear all;
close all;
clc;

tic
tiempo_a_graficar = 20;
save_interval = 1000;     % Guardar cada 1 unidad de tiempo (k=0.0001)

params.lambda_S = 0.01; params.lambda_E = 0.01; params.lambda_I = 0.01; 
params.lambda_Q = 0.01; params.lambda_R = 0.01; params.beta1  = 0.06;
params.beta2 = 0.07; params.mu = 0.06; params.delta = 0.05; 
params.gamma = 0.04; params.alfa = 0.05; params.rho = 0.03;
params.Lambda = 1; params.h = 0.01; params.k = 0.0001;

total_pasos = ceil(tiempo_a_graficar / params.k);
t_exp = 1:total_pasos;

%% Creación de la malla en t=0
S0 = ones(101, 101);  % CPU array para compatibilidad
E0 = zeros(101, 101);
I0 = zeros(101, 101);
Q0 = zeros(101, 101);
R0 = zeros(101, 101);

% Condiciones iniciales región central
S0(48:52,48:52) = 0.70;  % 70% susceptibles
E0(48:52,48:52) = 0.20;  % 20% expuestos
I0(48:52,48:52) = 0.10;  % 10% infectados

%% Simulación
[S_hist, E_hist, I_hist, Q_hist, R_hist] = modelo_SEIQR(t_exp, params, S0, E0, I0, Q0, R0, save_interval);

%% Graficar en 2D (tiempos específicos: 1, 20, 40, 60, 80, 100)
times_to_plot = [1, 20, 40, 60, 80, 100];
for t_val = times_to_plot
    figure('Name', ['Resultados 2D en t=', num2str(t_val)]);
    
    tiledlayout(3, 2, 'TileSpacing', 'compact');
    title(['Mapas de Calor en t=', num2str(t_val)]);
    colormap("jet");
    
    nexttile;
    imagesc(S_hist{t_val});
    set(gca, 'YDir', 'normal');
    clim([0, 1]);
    title('Susceptibles');
    colorbar;
    
    nexttile;
    imagesc(E_hist{t_val});
    set(gca, 'YDir', 'normal');
    clim([0, 1]);  % Ajustado para mejor visualización
    title('Expuestos');
    colorbar;
    
    nexttile;
    imagesc(I_hist{t_val});
    set(gca, 'YDir', 'normal');
    clim([0, 1]);  % Ajustado para mejor visualización
    title('Infectados');
    colorbar;
    
    nexttile;
    imagesc(Q_hist{t_val});
    set(gca, 'YDir', 'normal');
    clim([0, 1]);  % Ajustado para mejor visualización
    title('Cuarentenados');
    colorbar;
    
    nexttile;
    imagesc(R_hist{t_val});
    set(gca, 'YDir', 'normal');
    clim([0, 1]);
    title('Recuperados');
    colorbar;
end

t = toc;
fprintf('Tiempo total de ejecución: %.2f segundos\n', t);