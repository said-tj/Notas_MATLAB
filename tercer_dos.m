% ------------------------------------------------------------
% 1-D wave equation  u_tt = c^2 u_xx  (c=1, L=1)
% FEM lineal en x (3 elementos)  +  DIFERENCIA CENTRAL en t
% ------------------------------------------------------------
clc;  clear;

% ----- malla espacial
L  = 1;           % longitud
Ne = 3;           % nº de elementos ⇒ 4 nodos
h  = L/Ne;        % 1/3
x  = linspace(0,L,Ne+1).';

% ----- paso temporal
dt = 0.2;         % solicitado
T  = 2;           % tiempo final
Nt = round(T/dt); % nº de pasos

% ----- matrices ELEMENTO lineal
Me = h/6 * [2 1; 1 2];
Ke = 1/h * [ 1 -1 ; -1 1];     %   c = 1

% ----- ensamblaje global (4×4)
nNod = Ne+1;                % 4 nodos
M = zeros(nNod);
K = zeros(nNod);
for e = 1:Ne
    idx        = [e e+1];
    M(idx,idx) = M(idx,idx) + Me;
    K(idx,idx) = K(idx,idx) + Ke;
end

% ----- condición Dirichlet homogénea en nodos 1 y 4 (índices 1 y 4)
free = 2:3;                 % nodos interiores
Mr   = M(free,free);
Kr   = K(free,free);

A = dt^2 * (Mr\Kr);         %  matriz dinámica 2×2

% ----- condición inicial  u(x,0)=sin(pi x),  u_t(x,0)=0
U0      = sin(pi*x(2:3));   % [U1; U2] en t = 0
Um1     = U0;               % t = −dt  (porque velocidad inicial = 0)

Uh      = zeros(2,Nt+1);    % historial para graficar
Uh(:,1) = U0;

% ----- bucle en el tiempo (diferencia central)
Un   = U0;
Unm1 = Um1;
for n = 1:Nt
    Unp1 = 2*Un - Unm1 - A*Un;
    Unm1 = Un;
    Un   = Unp1;
    Uh(:,n+1) = Un;
end

% ----- reconstrucción incluyendo los ceros de contorno
X = x.';
Tvec = 0:dt:T;
Ufull = [ zeros(1,Nt+1) ; Uh ; zeros(1,Nt+1) ];  % 4×(Nt+1)

% ----- animación rápida
figure
for k = 1:length(Tvec)
    plot(X, Ufull(:,k), '-o','LineWidth',1.4);
    xlabel('x'), ylabel('u(x,t)')
    title(sprintf('t = %.1f s', Tvec(k)))
    axis([0 1 -1 1]), grid on, drawnow
end
