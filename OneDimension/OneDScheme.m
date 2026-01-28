function [S_hist, E_hist, I_hist, Q_hist, R_hist] = OneDScheme(t_exp,...
    params, S0, E0, I0, Q0, R0, save_interval)
    % Pre-alocación de celdas para guardar los resultados
    num_saves = floor(length(t_exp) / save_interval);
    S_hist = cell(1, num_saves);
    E_hist = cell(1, num_saves);
    I_hist = cell(1, num_saves);
    Q_hist = cell(1, num_saves);
    R_hist = cell(1, num_saves);
    
    % Inicialización de vectores de estado
    S_ant = S0; E_ant = E0; I_ant = I0; Q_ant = Q0; R_ant = R0;
    S_act = S0; E_act = E0; I_act = I0; Q_act = Q0; R_act = R0;
    
    phi = params.k;
    r = params.k / (params.h^2);
    save_counter = 1;

    for n = 1:length(t_exp)
        % --- CÁLCULO VECTORIZADO DE nu ---
        if n == 1
            nu_S = 0; nu_E = 0; nu_I = 0; nu_Q = 0; nu_R = 0;
        else
            nu_S = - (S_act - S_ant) ./ (params.k * (S_act + eps));
            nu_E = - (E_act - E_ant) ./ (params.k * (E_act + eps));
            nu_I = - (I_act - I_ant) ./ (params.k * (I_act + eps));
            nu_Q = - (Q_act - Q_ant) ./ (params.k * (Q_act + eps));
            nu_R = - (R_act - R_ant) ./ (params.k * (R_act + eps));
        end

        % --- ACTUALIZACIÓN VECTORIZADA DE TODAS LAS MATRICES ---
        S_num = S_act + params.lambda_S * r * neighborsOD(S_act) + phi * (params.Lambda + params.rho * Q_act);
        S_den = 1 + 2 * params.lambda_S * r * (1 + nu_S * params.k) + phi * (params.beta1 * E_act + params.beta2 * I_act + params.mu);
        S_sig = S_num ./ S_den;

        E_num = E_act + params.lambda_E * r * neighborsOD(E_act) + phi * params.beta1 .* E_act .* S_sig;
        E_den = 1 + 2 * params.lambda_E * r * (1 + nu_E * params.k) + phi * (params.delta + params.mu);
        E_sig = E_num ./ E_den;
        
        I_num = I_act + params.lambda_I * r * neighborsOD(I_act) + phi * (params.beta2 .* I_act .* S_sig + params.delta * E_sig);
        I_den = 1 + 2 * params.lambda_I * r * (1 + nu_I * params.k) + phi * (params.gamma + params.mu);
        I_sig = I_num ./ I_den;

        Q_num = Q_act + params.lambda_Q * r * neighborsOD(Q_act) + phi * params.gamma * I_sig;
        Q_den = 1 + 2 * params.lambda_Q * r * (1 + nu_Q * params.k) + phi * (params.alfa + params.rho + params.mu);
        Q_sig = Q_num ./ Q_den;

        R_num = R_act + params.lambda_R * r * neighborsOD(R_act) + phi * params.alfa * Q_sig;
        R_den = 1 + 2 * params.lambda_R * r * (1 + nu_R * params.k) + phi * params.mu;
        R_sig = R_num ./ R_den;
        
        % Actualización de estados para el siguiente paso
        S_ant = S_act; E_ant = E_act; I_ant = I_act; Q_ant = Q_act; R_ant = R_act;
        S_act = S_sig; E_act = E_sig; I_act = I_sig; Q_act = Q_sig; R_act = R_sig;

        % --- Guardado selectivo de datos ---
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