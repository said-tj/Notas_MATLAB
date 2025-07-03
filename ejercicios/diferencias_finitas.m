function diferencias_finitas()
    %{
    Ecuación:
        u_xx + u_yy = 50*sin(3*pi*x)*cos((pi/4)*y)
    en  0 <= x <= 1, 0 <= y <= 2

    Condiciones de frontera:
    1) u(0,y)   = 0                
    2) u(x,0)   = 0                
    3) du/dy=0  en x=1             
                                   
    4) du/dy=10*sin(pi*x) en y=2   
    
    Se usa un mallado de:
       n = 5 en x -> dx = 1/5
       m = 4 en y -> dy = 2/4

    %}

    %% 1. Parámetros de la malla
    nx = 5;       % número de subintervalos en x
    ny = 4;       % número de subintervalos en y
    dx = 1.0 / nx;  % 0.2
    dy = 2.0 / ny;  % 0.5

    x = 0 : dx : 1;   % 6 nodos en x:  i=0..5
    y = 0 : dy : 2;   % 5 nodos en y:  j=0..4

    % Cantidad de nodos interiores (i=1..4, j=1..3):
    Nx_int = 4;  % en x
    Ny_int = 3;  % en y
    N = Nx_int * Ny_int;  % total de incógnitas = 12

    %% Función f(x,y)
    f = @(xx,yy) 50 * sin(3*pi*xx) .* cos((pi/4)*yy);

    %% Nodos interiores
    idx = @(i,j) (j-1)*Nx_int + i; 

    % A y b para el sistema lineal
    A = zeros(N, N);
    b = zeros(N, 1);

    %% Ecuación en cada nodo interior
    % Las aproximación de 2a orden:
    %   (U_{i+1,j} - 2U_{i,j} + U_{i-1,j}) / dx^2 + 
    %   (U_{i,j+1} - 2U_{i,j} + U_{i,j-1}) / dy^2 = f(x_i, y_j).
    % i en [1..4], j en [1..3].

    for j = 1 : Ny_int   % j=1..3
        for i = 1 : Nx_int   % i=1..4

            % índice de fila en el sistema lineal
            row = idx(i,j);

            % coordenadas del nodo:
            Xi = x(i);   % i corre de 1..4 => x(i) 
            Yj = y(j);   % j corre de 1..3 => y(j)

            % El central: -2/dx^2 -2/dy^2
            % Empezamos sumando la parte para U_{i,j}
            A(row, row) = -2/(dx^2) - 2/(dy^2);

            %  U_{i+1,j} (hacia la derecha) --
            if (i+1 <= 4)
                % es nodo interior
                col = idx(i+1, j);
                A(row, col) = 1/(dx^2);
            else
                % i+1=5 => borde x=1, 
                % Usaremos la simplificación du/dy=0 => U(1,y)=0
                % => U_{5,j}=0 (Dirichlet 'forzado').
                % Entonces no suma nada a la matriz, sino a b si fuese !=0 
                % pero aquí es 0 => no hace nada.
            end

            %  U_{i-1,j} (hacia la izquierda) --
            if (i-1 >= 1)
                % es nodo interior
                col = idx(i-1, j);
                A(row, col) = 1/(dx^2);
            else
                % i-1=0 => borde x=0 => U_{0,j}=0
                % => no aporta término a A, ni a b (porque es 0)
            end

            % U_{i,j+1} (hacia arriba) --
            if (j+1 <= 3)
                % es nodo interior
                col = idx(i, j+1);
                A(row, col) = 1/(dy^2);
            else
                % j+1=4 => borde superior y=2 => Neumann: du/dy=10*sin(pi*x_i)
                % => U_{i,4} = U_{i,3} + (10 sin(pi*x_i))*dy = U_{i,3} + 5*sin(pi*x_i).
                %
                % => Efecto: U_{i,j+1} -> "U_{i,4} = U_{i,3} + 5*sin(pi*x_i)"
                % Sustitución: (U_{i,4} - 2U_{i,j} + U_{i,j-1})/dy^2 
                % se convierte en ( (U_{i,3}+5*sin(...)) - 2U_{i,j} + U_{i,j-1})/dy^2
                % => en la matriz, U_{i,4} se reemplaza por U_{i,3} => un 1/(dy^2) para U_{i,3}
                % y la parte "5*sin(...)/dy^2" va a b.
                col = idx(i, j);  % el propio U_{i,3}
                A(row, col) = A(row, col) + 1/(dy^2);

                % término extra en b:
                % + (5*sin(pi*x_i))/dy^2
                b(row) = b(row) + (5 * sin(pi*Xi)) / (dy^2);
            end

            %  U_{i,j-1} (hacia abajo) 
            if (j-1 >= 1)
                % es nodo interior
                col = idx(i, j-1);
                A(row, col) = A(row, col) + 1/(dy^2);
            else
                % j-1=0 => borde inferior => U_{i,0}=0
                % => no aporta nada a la matriz ni a b
            end

            % Término del lado derecho: f(x_i, y_j)
            b(row) = b(row) + f(Xi, Yj);
        end
    end

    %% Para resolver el sistema lineal
    Uvec = A \ b;   % (size 12 x 1)

    % Uvec(k) = U_{i,j}, con k=idx(i,j).

    U_interior = zeros(Nx_int, Ny_int);
    for j = 1 : Ny_int
        for i = 1 : Nx_int
            U_interior(i,j) = Uvec(idx(i,j));
        end
    end

    %% Reconstruir en toda la malla incluyendo los bordes
    % No olvidar que ppara graficar ordenadamente, creamos Ufull(0..5,0..4).
    Ufull = zeros(nx+1, ny+1);

    % Rellenamos nodos interiores
    for j = 1 : Ny_int    % j=1..3
        for i = 1 : Nx_int
            Ufull(i,j) = U_interior(i,j);
        end
    end

    % Borde x=0 => 0 
    % Borde y=0 => 0 
    % Borde x=5 => 0 
    % Borde y=4 => U_{i,4} = U_{i,3} + 5*sin(pi*x_i)
    for i = 1 : Nx_int
        X_i = x(i);
        Ufull(i,4) = Ufull(i,3) + 5 * sin(pi*X_i);
    end
    % La esquina x=5,y=4 => 0 (consistente con la simplificación adoptada).

    %% Graficar la solución
    [XX,YY] = meshgrid(x, y);  % meshgrid se invierte: X ~ columnas, Y ~ filas

    % Transponemos Ufull para que coincida con la convención de meshgrid:
    Uplot = Ufull';   % Uplot(j,i) = Ufull(i,j)

    figure
    surf(XX, YY, Uplot);
    xlabel('x'); ylabel('y'); zlabel('u(x,y)');
    title('Solución aproximada por diferencias finitas (mallado 5x4)');
end
