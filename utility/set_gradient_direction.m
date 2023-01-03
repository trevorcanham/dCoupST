function [negative_grad] = Set_gradient_direction(n,D,x0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the direction of the gradient in the normalization
%
% References: "NESTED NORMALIZATIONS FOR DECOUPLING GLOBAL FEATURES", 
% Javier Portilla & Eduardo Martinez-Enriquez. International Conference on
% Image Processing (ICIP) 2018. 
%
% INPUT:
%   n:                  Moment's order
%   D:                  Desired value
%   x0:                 x after normalization of the n-1 previous moments                        
%
% OUTPUT:
%   negative_grad:      Flag indicating the gradient direction
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


negative_grad=0;
Actual_feat=mean(x0.^n); 
if Actual_feat>D(n)
negative_grad=1;
end    


end

