% Limpiamos la consola y las variables.
clc;
clear all;

% Definimos una función.
function resultado = f_suma(a, b)
    resultado = a + b;
end

% Llamamos a la función y le damos argumentos.
respuesta = f_suma(2, 8);

% Mostramos al usuario el resultado.
%fprintf('La suma es %.2f\n',respuesta);

disp(['La suma es: ', num2str(resupuesta)]);
