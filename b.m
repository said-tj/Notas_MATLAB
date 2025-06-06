% =======================================================================
%  calor_FEM_Splines_Ecuaciones.m
% -----------------------------------------------------------------------
%  Ecuación de calor 1-D resuelta con *splines cúbicos uniformes* (FEM)
%  + θ-método en el tiempo.  El programa:
%     · Muestra fracciones de las integrales C_ij y K_ij.
%     · Imprime las matrices de masa (M) y rigidez (K).
%     · Lista TODAS las ecuaciones A·U^{n+1}=B·U^n en formato “cuaderno”.
%     · Evoluciona u(x,t) y grafica la solución.
%
%  Autor: Saïdt – Jun-2025
% =======================================================================

clear;  clc;  close all;
format short g;      % números compactos
format loose;

%% -------------------- DATOS DE ENTRADA --------------------------------
fprintf('=========== FEM (cubic splines) + θ-method ===========\n');
L      = input('Domain length  L                   = ');
alpha  = input('Diffusivity    α                   = ');
N      = input('Number of elements  N (≥4)         = ');
dt     = input('Time step      Δt                  = ');
Tfin   = input('Final time     T                   = ');
theta  = input('θ (1=implicit, 0.5=Crank-Nicolson) = ');

% Funciones inicial / frontera  (ajústalas si hace falta)
f  = @(x) sin(pi*x/L);    % u(x,0)
g0 = @(t) 0*t;            % u(0,t)
gL = @(t) 0*t;            % u(L,t)

%% -------------------- MALLA Y DOF -------------------------------------
h    = L/N;
x    = linspace(0,L,N+1)';      % nodos
nDOF = N-1;                     % incógnitas interiores U1…U_{N-1}

fprintf('\nSTEP 1  (Substitution)  —  DOF = %d\n',nDOF);
fprintf(  'STEP 2  (Galerkin weighting) — same basis\n');

%% ---------------- COEFICIENTES DEL PROFESOR ---------------------------
%   C_{|i-j|} · h     |     K_{|i-j|} · 1/h
c = [ 5/16 ,   3/112 , 129/2240 , 1/112 ];     % masa
k = [13/20 , -13/60 ,     1/60 ,     0 ];      % rigidez

fprintf('\nSTEP 3  (Integrals)\n');
fprintf('  |i-j|   C_ij (·h)         K_ij (/h)\n');
for m = 0:3
    fprintf('   %d      %s              %s\n',m,rat(c(m+1)),rat(k(m+1)));
end

% Diagonales para spdiags
mainC = c(1)*ones(nDOF,1); off1C=c(2)*ones(nDOF,1); off2C=c(3)*ones(nDOF,1); off3C=c(4)*ones(nDOF,1);
mainK = k(1)*ones(nDOF,1); off1K=k(2)*ones(nDOF,1); off2K=k(3)*ones(nDOF,1); off3K=k(4)*ones(nDOF,1);

M = h        * spdiags([off3C off2C off1C mainC off1C off2C off3C], -3:3, nDOF, nDOF);
K = (alpha/h)* spdiags([off3K off2K off1K mainK off1K off2K off3K], -3:3, nDOF, nDOF);

fprintf('\nMatrices interiores (banda 7)  —  size %d×%d\n',nDOF,nDOF);
if N<=12
    disp('M ='); disp(full(M));
    disp('K ='); disp(full(K));
else
    disp('   (N grande ⇒ impresión completa omitida)');
end

%% ---------------- CONDICIONES INICIALES -------------------------------
fprintf('\nSTEP 4  (IC & BC)\n');
U = f(x(2:end-1));               % vector inicial
disp('Initial U ='); disp(U.');

%% ---------------- SISTEMA LINEAL EN EL TIEMPO -------------------------
fprintf('\nSTEP 5  (System)   A·U^{n+1} = B·U^{n}\n');
A = M + theta    *dt*K;
B = M - (1-theta)*dt*K;

% Imprime las ecuaciones (convertimos A a denso para sprintf)
print_equations(full(A),'U_next');

% Factoriza para el bucle temporal
[Llu,Ulu] = lu(A);

%% ---------------- BUCLE TEMPORAL + GRÁFICA ----------------------------
nSteps = ceil(Tfin/dt);   t=0;
figure('Name','u(x,t) — FEM cubic splines');
for n = 1:nSteps
    t   = t + dt;
    rhs = B*U;
    U   = Ulu \ (Llu \ rhs);
    Uplot = [g0(t); U; gL(t)];
    plot(x,Uplot,'-o','LineWidth',1.4);
    xlabel('x'); ylabel('u(x,t)'); grid on;
    title(sprintf('time  t = %.4f',t)); ylim([-1 1]);
    drawnow;
end
fprintf('\nSimulation finished.\n');

%% ======================================================================
%%        AUXILIAR: Imprime ecuaciones fila-por-fila
%% ======================================================================
function print_equations(A,varName)
% A debe venir como matriz DENSO (full(A))
    n = size(A,1);
    fprintf('\nPrinting %d equations (coefficients 6-sig-fig):\n\n',n);
    for i = 1:n
        row = A(i,:);           % row es vector double
        eq  = sprintf('Eq %-3d: ',i);
        first = true;
        for j = 1:n
            a = row(j);
            if abs(a) < 1e-12,  continue,  end
            if ~first
                eq = [eq, sprintf(' %c ', char('+' *(a>0) + '-'*(a<0)))];
            elseif a < 0
                eq = [eq, '-'];  % signo inicial negativo
            end
            a = abs(a);
            if abs(a-1) > 1e-12
                eq = [eq, sprintf('%.6g*%s%d',a,varName,j)];
            else
                eq = [eq, sprintf('%s%d',varName,j)];
            end
            first = false;
        end
        eq = [eq, sprintf(' = rhs%d',i)];
        disp(eq);
    end
    fprintf('\n');
end
