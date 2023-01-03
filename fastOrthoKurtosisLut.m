function [ok,xout,t0,lut] = fastOrthoKurtosisLut(x,lut0)
% [ok,xout,t0] = fastOrthoKurtosis(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    This function obtains the ortho-kurtosis from a quasi-analytic
%    solution. See Section 2.3 of the paper for details.
%
%    References: "DETERMINISTIC FEATURE DECOUPLING BY SURFING INVARIANCE MANIFOLDS", 
%    Eduardo Martinez-Enriquez & Javier Portilla. International Conference on Acoustics, 
%    Speech, and Signal Processing (ICASSP) 2020. 
%
% Authors:
%  - Javier Portilla & Eduardo Martinez-Enriquez <emenriquez@tsc.uc3m.es>
% 
%     Copyright (C)  2020 Javier Portilla & Eduardo Martinez-Enriquez .
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

[Ny,Nx] = size(x);
x = x(:) - mean(x(:));

mu3 = mu3_central(0,x);
if mu3>0
    t_guess = [(1-1e-3)/min(x) 0];
else    
    t_guess = [0 (1-1e-3)/max(x)];
end

t0 = fzero(@(t) mu3_central(t,real(x)),t_guess);


%Impose 0 skew
x1 = x./(1-t0*x);
lut1 = lut0./(1-t0*lut0);

% Measure (ortho)kurtosis on the 0-skew sample
x2 = x1 - mean(x1);
lut2 = lut1 - mean(x1); 
x22 = x2.^2;
sigma2 = mean(x22);
x22 = x22/sigma2;
ok = mean(x22.^2);

if nargout>1
    xout = x2.'/sqrt(sigma2); % normalized output: zero-skew, unit-variance, zero-mean
    lut = lut2./sqrt(sigma2);
    xout = reshape(xout, Ny, Nx);
end    

