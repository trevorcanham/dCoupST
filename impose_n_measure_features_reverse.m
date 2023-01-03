function [measured_central_moments,imposed_features,x] = impose_n_measure_features_reverse(x,desired_features, rel_tol)

%  [measured_central_moments,measured_new_features,imposed_features,x] = impose_n_measure_features_reverse(x,desired_features, rel_tol)
%
% This is a modified version of impose_n_measure_features that takes a signal that
% has been progressively normalized, and it does the exact reverse transformation,
% for imposing a new desired set of values for the orthomoments.
% See also, besides impose_n_measure_features.m,
% impose_n_measure_features_synth.m (the former only does the forward normalization, 
% whereas the latter does both forward and reverse normalization).
% See, besides the reference below (IEEE ICIP 2018) also IEEE MMSP 2020 ("Controlled
% Feature Adjustment...") by the same authors.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script starts with a given set of data defined in a vector x and, following 
% Portilla and Martinez's method, sequentially changes x until it has a sample mean, variance, 
% skewness and kurtosis defined by the user.
% The function uses matlab ODE solvers that adaptively choses the integration step size.
%
% References: "NESTED NORMALIZATIONS FOR DECOUPLING GLOBAL FEATURES", 
% Javier Portilla & Eduardo Martinez-Enriquez. International Conference on
% Image Processing (ICIP) 2018. 
%
% INPUT:
%   x:                  Signal to be measured/processed
%   desired_features:   Vector with the desired values of the moments after
%                       normalization (by default: Gaussian)
%                       (alternatively: Number N of moments, for Gaussian
%                       reference).
%  rel_tol:             Relative tolerance for the solution of the ODE and thus 
%                       for the adjustment of the parameters' values
%
% OUTPUT:
%   measured_central_moments:   Moments of the signal before normalization
%   measured_new_features:      Moments of the signal after progressive normalizations
%   x:                          Transformed signal after nested normalizations


%
% Authors:
%  - Javier Portilla & Eduardo Martinez-Enriquez <emenriquez@tsc.uc3m.es>. 
%    Marzo 2018, IO-CSIC, Madrid.
% 
%     Copyright (C)  2018 Javier Portilla & Eduardo Martinez-Enriquez .
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if exist('desired_features','var'),
    N = length(desired_features);
    if N==1,
        N = desired_features;
    end    
else
    warning('Number of measured features assumed to be 8 (N+1)')
    N = 7;
end    
[Ny,Nx] = size(x);
Nel = numel(x);
x = reshape(x,Nel,1);

% Measure the standardized central moments of the signal before
% transformation

media = mean(x);
mcm(1) = media;
variance = mean(x.^2) - media^2;
mcm(2) = variance;
if N>2,
%for n = 3:N+1,
for n = 3:N
    mcm(n) = mean(((x-media)/sqrt(variance)).^n);
end    
end


measured_central_moments = mcm;

% Compute reference vector D for the moments of a univariate zero-mean
% Gaussian density.

% Normalization values
    a = 1;
    for n = 1:N,
        if n/2 ~= round(n/2),
           Dn(n) = 0; % Even symmetry, odd moments are all zero.
        else
           a = (n-1)*a; % 2->1,4->3 (3x1),6->15 (5x3x1), 8->105 (7x5x3x1), 10->945 (9x7x5x3x1)....
           Dn(n) = a;
        end
    end
    
SNR_th = 150; % Threshold in dB in SNR w.r.t. the moments of a normalized signal, to consider that the signal is acceptable    
snrr = snr(mcm(1:N-1),mcm(1:N-1)-Dn(1:N-1));
if snrr<SNR_th
    error('The passed signal has not been properly normalized.')
    return
end    
    
% Desired features
D = desired_features;

% Set tolerance values for the ODE solution
if ~exist('rel_tol'),
    rel_tol = 1e-7;
end    

x0 = x;

% the orthofeatures are imposed in reverse order

if N>2
for n = N:-1:3

[negative_grad] = set_gradient_direction(n,D,x0);

warning off
% Solve the ODE
%NN = 10000; % I do not see any impact of this parameter....
%NN = 1000;
NN = 100;
%NN = 10;
%NN = 1;
opt = odeset('Events', @(t,x) myEvent_n(t,x,n,negative_grad,D),'RelTol',rel_tol,'AbsTol',rel_tol*1e-4);
opt.MaxStep
[t,x] = ode45(@(t,x) grad_n(t,x,n,negative_grad),[0:NN],x0, opt);
%size(x)
x0 = x(end,:);

gn = mean(x0.^n); 
mnf(n) = gn;

end % for    
end % if N>2

% The first two normalizations are done analytically in one shot

media = mean(x0);
variance = mean((x0-media).^2);
x0 = (x0-media)*sqrt(D(2)/variance) + D(1);
media = mean(x0);
variance = mean((x0-media).^2);

mnf(1) = media;
mnf(2) = variance;

imposed_features = mnf;


x = reshape(x0,Ny,Nx);







