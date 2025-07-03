% ============================================================
%    calor_FEM_SplinesCubicOSalidaCompleta.m
% ------------------------------------------------------------
%  FEM 1-D con *splines cúbicos básicos uniformes*  +  θ-método.
%  Imprime:
%    • Integrales C_ij y K_ij (forma exacta y numérica)
%    • Matrices M y K completas y recortadas
%    • Cada uno de los 5 pasos pedidos
%  Autor: Saïdt · Jun-2025
% ============================================================

clear;  clc;  close all;
format rational   % para que MATLAB muestre fracciones cuando pueda

%% ======  Datos del problema ======
fprintf('=== Difusión 1-D con splines cúbicos — Salida completa ===\n');
L      = input('Longitud del dominio L              = ');
alpha  = input('Difusividad      alpha              = ');
N      = input('Número de elementos N (p.ej. 8–40)  = ');
dt     = input('Paso temporal     dt                = ');
Tfin   = input('Tiempo final      T                 = ');
theta  = input('θ-método  (1=implícito, 0.5=C-N)    = ');

% Funciones de datos (ajusta a tu caso)
f  = @(x) sin(pi*x/L);     % u(x,0)
g0 = @(t) 0*t;             % Dirichlet en x=0
gL = @(t) 0*t;             % Dirichlet en x=L

%% ======  Definición de malla y DOF ======
h   = L/N;                       % paso
x   = linspace(0,L,N+1)';        % nodos globales
nDOF= N-1;                       % nodos interiores (1…N-1)

fprintf('\nPaso 1) Sustitución  u ≈ w = Σ U_i(t) S_i(x)\n');
fprintf('        -> nDOF = %d grados de libertad (sin fronteras)\n',nDOF);

%% ======  Paso 2) Ponderación (Galerkin)  ======
fprintf('\nPaso 2) Ponderamos con las mismas funciones S_k(x) (Galerkin)\n');

%% ======  Paso 3) Integración: Construcción de C_ij y K_ij ======
fprintf('\nPaso 3) Integrales sobre el dominio — Matrices de masa y rigidez\n');

% ----------  COEFICIENTES DEL PROFESOR  ----------------------
%   C_{|i-j|}  multiplican  h
c0 =  5/16 ;    c1 =  3/112 ;   c2 = 129/2240 ;   c3 = 1/112 ;
%   K_{|i-j|}  multiplican  1/h
k0 =  13/20 ;   k1 = -13/60 ;   k2 =   1/60 ;    k3 = 0 ;

fprintf('\n  >> Integrales exactas (fracción · h  ó  fracción / h)\n');
fprintf('     |i-j|   C_(ij)               K_(ij)\n');
fprintf('       0     %s·h        %s/h\n',rat(c0),rat(k0));
fprintf('       1     %s·h        %s/h\n',rat(c1),rat(k1));
fprintf('       2     %s·h        %s/h\n',rat(c2),rat(k2));
fprintf('       3     %s·h        %s/h\n\n',rat(c3),rat(k3));

%  → vectores diagonales para spdiags
mainC = c0*ones(nDOF,1);  off1C=c1*ones(nDOF,1); off2C=c2*ones(nDOF,1); off3C=c3*ones(nDOF,1);
mainK = k0*ones(nDOF,1);  off1K=k1*ones(nDOF,1); off2K=k2*ones(nDOF,1); off3K=k3*ones(nDOF,1);

% Matrices globales interiores (ya sin nodos Dirichlet)
M = h * spdiags([off3C off2C off1C mainC off1C off2C off3C], -3:3, nDOF, nDOF);
K = (alpha/h) * spdiags([off3K off2K off1K mainK off1K off2K off3K], -3:3, nDOF, nDOF);

% ------  Mostrar matrices (formato denso para N≲12) ----------
fprintf('  Matriz de masa M (interior) — tamaño %d × %d\n',size(M));
if N<=12, disp(full(M)); else disp('  (matriz muy grande, se omite impresión completa)'); end
fprintf('\n  Matriz de rigidez K (interior) — tamaño %d × %d\n',size(K));
if N<=12, disp(full(K)); else disp('  (matriz muy grande, se omite impresión completa)'); end

%% ======  Paso 4) Condiciones iniciales y frontera ======
fprintf('\nPaso 4) Impongo condiciones iniciales y Dirichlet homogéneas\n');
U = f(x(2:end-1));             % vector inicial
disp('  Vector U(t=0):'); disp(U.');

%% ======  Paso 5) Sistema lineal en el tiempo ================
fprintf('\nPaso 5) θ-método → (M+θΔtK)U^{n+1} = (M-(1-θ)ΔtK)U^n\n');

A = M + theta*dt*K;            % lado izq.
B = M - (1-theta)*dt*K;        % lado der.
[Llu,Ulu] = lu(A);             % factoriza una sola vez
fprintf('  Factorización LU completa.\n');

%% --------   Bucle temporal con gráfica en vivo  -------------
nSteps = ceil(Tfin/dt); t=0;
figure('Name','Evolución u(x,t) — FEM splines cúbicos');
for n=1:nSteps
    t = t + dt;
    rhs = B*U;
    U   = Ulu \ (Llu \ rhs);
    % reconstruir con Dirichlet para graficar
    Uplot = [g0(t); U; gL(t)];
    plot(x,Uplot,'-o','LineWidth',1.4);
    xlabel('x'); ylabel('u(x,t)'); grid on;
    title(sprintf('Solución FEM — t = %.4f',t)); ylim([-1 1]); drawnow;
end
fprintf('\n>>> Simulación terminada sin errores.\n');

% ====================  FIN DEL SCRIPT ========================
