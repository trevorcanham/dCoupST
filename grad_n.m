
function [dydt] = grad_n(t,x,n,negative_grad)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE system to be integrated. Analytical expressions of the gradient of
% order n
%
% References: "NESTED NORMALIZATIONS FOR DECOUPLING GLOBAL FEATURES", 
% Javier Portilla & Eduardo Martinez-Enriquez. International Conference on
% Image Processing (ICIP) 2018. 
%
% INPUT:
%   n:                  Moment's order
%   D:                  Desired value
%   x:                  x after normalization of the n-1 previous moments 
%   negative_grad:      Flag indicating the gradient direction
%
% OUTPUT:
%   dydt:               Expressions of the gradients to be integrated                   
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





mu = zeros(2*(n-1)-1,1);
N = length(x);
ac = ones(length(x),1);
xp = zeros(length(x),max(2*(n-1)-1,n));
for j = 1:max(2*(n-1)-1,n),
    ac = ac.*x;
    xp(:,j) = ac;
    mu(j) = mean(ac);
end    

if n>2,
for i = 1:n-2
    for j = i:n-1
        c(1,i,j) = mu(i+j) - mu(i)*mu(j);
    end
end    
end

if n>3,
for k = 2:n-2,
for i = 2:n-2
    for j = i:n-1
        c(k,i,j) = c(k-1,i,j)*c(k-1,k-1,k-1) - c(k-1,k-1,i)*c(k-1,k-1,j);
    end
end    
end
end

k = 1;
for j = k+1:n,
    a(j,k) = mu(j-1);
end

if n>2,
for k = 2:n-1,
    for j = k+1:n,
        a(j,k) = c(k-1,k-1,j-1)/c(k-1,k-1,k-1);
    end
end    
end

g = zeros(length(x),n);
g(:,1) = ones(length(x),1);
for i = 2:n,
    g(:,i) = xp(:,i-1);
    for j = 1:i-1,
        g(:,i) = g(:,i) - a(i,j)*g(:,j);
    end
end    
    

% g1 = 1;
% g2 = x - a21*g1;
% g3 = x.^2 - a31*g1 - a32*g2;
% g4 = x.^3 - a41*g1 - a42*g2 - a43*g3;
% g5 = x.^4 - a51*g1 - a52*g2 - a53*g3 - a54*g4;

dydt = g(:,n);
if negative_grad==0
% dydt=dydt/norm(dydt);
dydt;
else
% dydt=-dydt/norm(dydt);   
dydt=-dydt;
end
% gn = g(:,n);
% mun = mu(n);
