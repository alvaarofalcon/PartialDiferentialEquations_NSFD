function [S_hist, E_hist, I_hist, Q_hist, R_hist] = modelo_SEIQR(t_exp, params, S0, E0, I0, Q0, R0, save_interval)
    % Pre-alocación de celdas
    num_saves = floor(length(t_exp) / save_interval);
    S_hist = cell(1, num_saves);
    E_hist = cell(1, num_saves);
    I_hist = cell(1, num_saves);
    Q_hist = cell(1, num_saves);
    R_hist = cell(1, num_saves);
    
    % Inicialización
    S_ant = S0; E_ant = E0; I_ant = I0; Q_ant = Q0; R_ant = R0;
    S_act = S0; E_act = E0; I_act = I0; Q_act = Q0; R_act = R0;
    
    phi = params.k;
    r = params.k / (params.h^2);
    save_counter = 1;
    
    % Kernel para vecinos (difusión)
    kernel = [0 1 0; 1 0 1; 0 1 0];  % Vecinos: norte, sur, este, oeste

    for n = 1:length(t_exp)
        % --- Cálculo de ν (para n>=2) ---
        if n == 1
            nu_S = 0; nu_E = 0; nu_I = 0; nu_Q = 0; nu_R = 0;
        else
            nu_S = - (S_act - S_ant) ./ (params.k * S_act + 1e-9);
            nu_E = - (E_act - E_ant) ./ (params.k * E_act + 1e-9);
            nu_I = - (I_act - I_ant) ./ (params.k * I_act + 1e-9);
            nu_Q = - (Q_act - Q_ant) ./ (params.k * Q_act + 1e-9);
            nu_R = - (R_act - R_ant) ./ (params.k * R_act + 1e-9);
        end
        
        % --- Cálculo de vecinos (difusión) ---
        neighbors_S = conv2(S_act, kernel, 'same');
        neighbors_E = conv2(E_act, kernel, 'same');
        neighbors_I = conv2(I_act, kernel, 'same');
        neighbors_Q = conv2(Q_act, kernel, 'same');
        neighbors_R = conv2(R_act, kernel, 'same');
        
        % --- Actualización SECUENCIAL (ecuación 9 del paper) ---
        % 1. Actualizar SUSCEPTIBLES (S)
        S_num = S_act + params.lambda_S * r * neighbors_S + phi * (params.Lambda + params.rho * Q_act);
        S_den = 1 + 4 * params.lambda_S * r * (1 + nu_S * params.k) + phi * (params.beta1 * E_act + params.beta2 * I_act + params.mu);
        S_act = S_num ./ S_den;
        S_act = aplicar_neumann(S_act);  % Condiciones de contorno
        
        % 2. Actualizar EXPUESTOS (E) - usa S actualizado
        E_num = E_act + params.lambda_E * r * neighbors_E + phi * params.beta1 .* E_act .* S_act;
        E_den = 1 + 4 * params.lambda_E * r * (1 + nu_E * params.k) + phi * (params.delta + params.mu);
        E_act = E_num ./ E_den;
        E_act = aplicar_neumann(E_act);
        
        % 3. Actualizar INFECTADOS (I) - usa S y E actualizados
        I_num = I_act + params.lambda_I * r * neighbors_I + phi * (params.beta2 .* I_act .* S_act + params.delta * E_act);
        I_den = 1 + 4 * params.lambda_I * r * (1 + nu_I * params.k) + phi * (params.gamma + params.mu);
        I_act = I_num ./ I_den;
        I_act = aplicar_neumann(I_act);
        
        % 4. Actualizar CUARENTENADOS (Q) - usa I actualizado
        Q_num = Q_act + params.lambda_Q * r * neighbors_Q + phi * params.gamma * I_act;
        Q_den = 1 + 4 * params.lambda_Q * r * (1 + nu_Q * params.k) + phi * (params.alfa + params.rho + params.mu);
        Q_act = Q_num ./ Q_den;
        Q_act = aplicar_neumann(Q_act);
        
        % 5. Actualizar RECUPERADOS (R) - usa Q actualizado
        R_num = R_act + params.lambda_R * r * neighbors_R + phi * params.alfa * Q_act;
        R_den = 1 + 4 * params.lambda_R * r * (1 + nu_R * params.k) + phi * params.mu;
        R_act = R_num ./ R_den;
        R_act = aplicar_neumann(R_act);
        
        % --- Guardar estados anteriores para siguiente paso ---
        S_ant = S_act;
        E_ant = E_act;
        I_ant = I_act;
        Q_ant = Q_act;
        R_ant = R_act;
        
        % --- Guardar resultados ---
        if mod(n, save_interval) == 0
            S_hist{save_counter} = S_act;
            E_hist{save_counter} = E_act;
            I_hist{save_counter} = I_act;
            Q_hist{save_counter} = Q_act;
            R_hist{save_counter} = R_act;
            save_counter = save_counter + 1;
        end
    end
end

function u = aplicar_neumann(u)
    % Condiciones de contorno Neumann
    u(1, :) = u(2, :);     % Borde superior
    u(end, :) = u(end-1, :); % Borde inferior
    u(:, 1) = u(:, 2);      % Borde izquierdo
    u(:, end) = u(:, end-1); % Borde derecho
end