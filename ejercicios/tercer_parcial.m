% ------------------------------------------------------------
% 1-D wave equation (c = 1, L = 1) with FEM (3 elems) + central time
% ------------------------------------------------------------
clear; clc;

% ---- parámetros
L  = 1;
c  = 1;
Ne = 3;           % número de elementos → 4 nodos
h  = L/Ne;        % 1/3
dt = 0.2;         % paso temporal
T  = 2;           % tiempo total a simular
Nt = round(T/dt); % nº de pasos

% ---- matrices de elemento (lineales)
Me = h/6 * [2 1; 1 2];
Ke = c^2/h * [1 -1; -1 1];

% ---- ensamblaje global (4 × 4)
M = zeros(4);
K = zeros(4);
for e = 1:Ne
    idx = [e e+1];
    M(idx,idx) = M(idx,idx) + Me;
    K(idx,idx) = K(idx,idx) + Ke;
end

% ---- reducir por Dirichlet (nodos 0 y 3 fijados)
free = [2 3];         % índices Matlab (nodos 1 y 2 físicos)
Mr = M(free,free);
Kr = K(free,free);

% ---- matriz A = dt^2 Mr^{-1} Kr
A = (dt^2) * (Mr\Kr);

% ---- condiciones iniciales
U0 = [ sin(pi/3) ; sin(2*pi/3) ];   % desplazamiento en t = 0
Um1 = U0;                           % t = -dt  (velocidad = 0)

% ---- pre-asignación para guardado
Uhist = zeros(2,Nt+1);
Uhist(:,1) = U0;

% ---- integración en el tiempo
Un = U0;
Unm1 = Um1;
for n = 1:Nt
    Unp1 = 2*Un - Unm1 - A*Un;   % esquema (∗)
    % avanzar
    Unm1 = Un;
    Un   = Unp1;
    Uhist(:,n+1) = Un;
end

% ---- reconstruir las 4 coordenadas (añadiendo los ceros de frontera)
tvec = 0:dt:T;
xnod = [0, 1/3, 2/3, 1];

figure
for k = 1:length(tvec)
    plot(xnod, [0; Uhist(:,k); 0], '-o');
    xlabel('x'), ylabel('u(x,t)')
    title(sprintf('t = %.2f s', tvec(k)))
    axis([0 1 -1 1]), grid on
    drawnow
end
