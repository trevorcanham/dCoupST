function im = perturbate_for_moments(im0, Q);

%im = perturbate_for_moments(im0, Q);
% Q: number of quantization levels

% Función para generar una perturbación entre -1/2 y 1/2-1/(2N)
% en un array 2D.
% La idea es que la perturbación sea discreta y única en cada píxel, 
% por lo que asegure que los valores resultantes sean todos diferentes.
% Además la perturbación es "discreta" también en el sentido de tener un
% bajo impacto perceptual, porque es una rampa, con lo que no genera saltos
% locales (es lo más suave que puede ser, en el sentido de,
% aproximadamente, minimizar la norma L1 de las diferencias de un pixel con
% sus vecinos, sujeto a la restricción de dar N valores diferentes 
% uniformemente equiespaciados).
% Además, la perturbación se escala adecuadamente para ser lo mayor posible
% sin afectar al resultado de su cuantización con los Q niveles originales
% (e.g., Q = 2^B, B = 8 bits), es decir, im0 = q(im) (perturbación
% reversible vía quantificación), por lo que no hay pérdida de información.
%
% El algortimo se basa en generar una rampa lineal que no repita valores y
% luego asignar valores de perturbación correspondientes al orden de mayor
% a menor. Esto resulta en una rampa sigmoidal con valores enteros no
% repetidos.
%
% Esta función está motivada en el contexto de ajustes con caracteristicas
% ortogonales ("orthofeatures"), para evitar que la multiplicidad de
% valores repetidos generen singularidades (puntos de ensillamiento, o
% "sillas") que hacen que el ajuste en las EDOs "se atasque", por generar dominios 
% de atracción a soluciones estacionarias espúreas (a valores de rango
% intermedios).
%
% Javier Portilla, Instituto de Óptica, CSIC, Abril 2021


% genera un array
[Ny,Nx] = size(im0);
[x,y] = meshgrid(1:Ny,1:Nx);

rep = true;

while rep 
%Primero genera un valor pseudo-random entre 0 y 1
r = rand;

% y el plano inclinado
v = x + y*r;

%close all; 
%figure
%imagesc(v); axis('image'); colorbar

v0 = v;
[sorted ind] = sort(v(:));

rep = (length(sorted)~=Ny*Nx); % comprueba que todos los píxeles son diferentes
end % while

ranking = zeros(1,Ny*Nx);
ranking(ind) = 1:Ny*Nx; % Obtain the ranking of each pixel
ranking = reshape(ranking,Ny,Nx);

%figure
%imagesc(ranking); axis('image'); colorbar

%Q = 2^8; % 8 bits de quantificación
reduction_factor = 0.95; % for preventing the perturbation to get out of the quantization interval
pert = 1/Q*(ranking-1-Ny*Nx/2)/(Ny*Nx); % mapea el ranking al rango reduction_factor*(-1/2,1/2)/Q
pert = reduction_factor*pert;

im = im0 + pert;