% wave_fem_3elements.m
% 1-D wave equation, c = 1, L = 1, Dirichlet homogéneo
% malla: 3 elementos lineales (nodos 0..3)
clear;  close all;  clc;

%% 0. Parámetros del problema
L      = 1;         % longitud
c      = 1;         % velocidad
Nelem  = 3;         % elementos (P1)
dt     = 0.2;       % paso de tiempo
Tfin   = 1.0;       % tiempo final (puedes ampliarlo)
nx     = Nelem+1;   % nodos totales  (0..3)
Nint   = nx-2;      % grados de libertad internos (i=1,2)
x      = linspace(0,L,nx);

%% 1. Matrices M y K (P1, malla uniforme h = L/3)
h   = L/Nelem;
M   = (h/6)*[4 1; 1 4];
K   = (1/h)*[ 1 -1; -1 1];
Minv = inv(M);                          % inversa pequeña 2×2
A    = 2*eye(2) - (dt^2)*Minv*K;        % matriz de avance explícito

%% 2. Condiciones iniciales u(x,0)=sin(pi x),  ut(x,0)=0
U0 = [sin(pi*x(2)); sin(pi*x(3))];      % valores en nodos internos
U1 = U0;                                % primer paso: velocidad cero  →  U1 = U0
Ut0 = zeros(2,1);                       % (no se usa en el esquema)

%% 3. Almacenamos resultados
nSteps  = round(Tfin/dt);
Unum    = zeros(2,nSteps+1);            % 2 DOF × (nSteps+1)
Unum(:,1) = U0;
Unum(:,2) = U1;

%% 4. Bucle temporal (diferencias centradas explícitas)
for n = 2:nSteps
    Unum(:,n+1) = A*Unum(:,n) - Unum(:,n-1);
end

%% 5. Añadimos los nodos de contorno (u = 0) para graficar
UnumFull = [zeros(1,nSteps+1);          % nodo 0
            Unum;                      % nodos 1 y 2
            zeros(1,nSteps+1)];        % nodo 3

%% 6. Solución analítica en la malla
t      = 0:dt:Tfin;
Uexact = ( sin(pi*x.') * cos(pi*t) );   % matriz (nodos × tiempos)

%% 7. Animación / gráficos
figure;  tiledlayout(1,2,"TileSpacing","compact")

% --- (a) Evolución numérica ---
nexttile;
for n = 1:nSteps+1
    plot(x,UnumFull(:,n),'b-o','LineWidth',1.4);  hold on
    plot(x,Uexact(:,n),'r--','LineWidth',1.2);
    hold off
    ylim([-1 1]); grid on
    title(sprintf('t = %.1f s',t(n)));
    legend('FEM (3 elems)','Exacta','Location','southoutside');
    xlabel('x'); ylabel('u(x,t)');
    drawnow;
end
% --- (b) Error en norma L2 vs. tiempo ---
err = sqrt( h*sum( (UnumFull-Uexact).^2 ,1) ); % ∫|e|² dx ≈ h Σ nodos
nexttile;
plot(t,err,'k-o','LineWidth',1.4);
xlabel('t'); ylabel('\|error\|_{L2}');
title('Error L2 frente a la exacta');
grid on
