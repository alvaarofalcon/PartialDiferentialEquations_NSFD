function Snbr = neighborsOD(S)
% Snbr = suma_primeros_vecinos_neumann(S)
% Devuelve para cada componente j la suma de sus primeros vecinos:
% - interior: S_j-1 + S_j+1
% - borde izquierdo: 2*S(2)  (Neumann: replica espejo)
% - borde derecho:  2*S(n-1)

n = numel(S);

% construir vector
Snbr = zeros(size(S));
% borde izquierdo
Snbr(1) = 2 * S(2);
% interior
Snbr(2:n-1) = S(1:n-2) + S(3:n);
% borde derecho
Snbr(n) = 2 * S(n-1);

end
