function im = perturbate_for_moments(im0, Q);

%im = perturbate_for_moments(im0, Q);
% Q: number of quantization levels

% Funci�n para generar una perturbaci�n entre -1/2 y 1/2-1/(2N)
% en un array 2D.
% La idea es que la perturbaci�n sea discreta y �nica en cada p�xel, 
% por lo que asegure que los valores resultantes sean todos diferentes.
% Adem�s la perturbaci�n es "discreta" tambi�n en el sentido de tener un
% bajo impacto perceptual, porque es una rampa, con lo que no genera saltos
% locales (es lo m�s suave que puede ser, en el sentido de,
% aproximadamente, minimizar la norma L1 de las diferencias de un pixel con
% sus vecinos, sujeto a la restricci�n de dar N valores diferentes 
% uniformemente equiespaciados).
% Adem�s, la perturbaci�n se escala adecuadamente para ser lo mayor posible
% sin afectar al resultado de su cuantizaci�n con los Q niveles originales
% (e.g., Q = 2^B, B = 8 bits), es decir, im0 = q(im) (perturbaci�n
% reversible v�a quantificaci�n), por lo que no hay p�rdida de informaci�n.
%
% El algortimo se basa en generar una rampa lineal que no repita valores y
% luego asignar valores de perturbaci�n correspondientes al orden de mayor
% a menor. Esto resulta en una rampa sigmoidal con valores enteros no
% repetidos.
%
% Esta funci�n est� motivada en el contexto de ajustes con caracteristicas
% ortogonales ("orthofeatures"), para evitar que la multiplicidad de
% valores repetidos generen singularidades (puntos de ensillamiento, o
% "sillas") que hacen que el ajuste en las EDOs "se atasque", por generar dominios 
% de atracci�n a soluciones estacionarias esp�reas (a valores de rango
% intermedios).
%
% Javier Portilla, Instituto de �ptica, CSIC, Abril 2021


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

rep = (length(sorted)~=Ny*Nx); % comprueba que todos los p�xeles son diferentes
end % while

ranking = zeros(1,Ny*Nx);
ranking(ind) = 1:Ny*Nx; % Obtain the ranking of each pixel
ranking = reshape(ranking,Ny,Nx);

%figure
%imagesc(ranking); axis('image'); colorbar

%Q = 2^8; % 8 bits de quantificaci�n
reduction_factor = 0.95; % for preventing the perturbation to get out of the quantization interval
pert = 1/Q*(ranking-1-Ny*Nx/2)/(Ny*Nx); % mapea el ranking al rango reduction_factor*(-1/2,1/2)/Q
pert = reduction_factor*pert;

im = im0 + pert;