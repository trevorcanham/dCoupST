function [w] = impose_orthofilters_BP2(hat_x,orthofeat,LHs,LLHs,LLLHs,LLLLs)
% This function imposes a set of orthofeature values to a previously
% normalized (up to the level N-1) image.
% The features, in their hierarchical order, are:
% 1) the mean
% 2) MSV (mean Square value) in the pixel domain
% 3) MSV at the output of a high-frequency band-pass filter, LH
% 4) MSV at the output of a medium-frequency band-pass filter, LLH
% 5) MSV at the output of a low-frequency band-pass filter, LLLH
% 6) MSV at the output of a very low-pass filter LLLL
%
% INPUT:
%   hat_x: normalized x, up to the N-1 (third) feature.
%   orthofeat: desired set of orthofeatures
%   (IMPLICIT, through global variables:)
%   LH, LLH, LLLH, LLLL:  spectral (centered) response of filters
%   LH_ref, LLH_ref, LLLH_ref 
%
% OUTPUT:
%   w:  transformed image (gray level)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Yf
global LH_ref
global LLH_ref
global LLLH_ref

global LH
global LLH
global LLLH
global LLLL

global LH_desired
global LLH_desired
global LLLH_desired
global LLLL_desired

LH = LHs;
LLH = LLHs;
LLLH = LLLHs;
LLLL = LLLLs;


% Now we go in reverse order, imposing the desired values (measured in
% another image or fabricated)

% First we correct the MSV output of the LL filter in hat_x, keeping
% the previous features with their references values:

[Ny,Nx] = size(hat_x);
if (2*round(Ny/2)~=Ny)||(2*round(Nx/2)~=Nx)
    error('Input image must have even size dimensions');
    return;
end


Yf = fftshift(fft2(hat_x));
LLLL_desired = orthofeat(6);

t0 = [0 0 0 0]; %initial guess
fun2 = @compute_f3_n_f4_n_f5_n_f6_of_normalized_x_BP;
options = optimoptions('fsolve','FunctionTolerance',1.0e-10,'OptimalityTolerance',1e-8);
%options = optimoptions('fsolve','FunctionTolerance',1.0e-6,'OptimalityTolerance',1e-6);
hat_t = fsolve(fun2,t0,options);

% We transform accordingly x
% First, we test that hat_f3 has been normalized:
Yft = Yf.* exp((abs(LH).^2*hat_t(1) + abs(LLH).^2*hat_t(2) + abs(LLLH).^2*hat_t(3) + abs(LLLL).^2*hat_t(4)));
Yft_s = Yft/sqrt(mean(abs(Yft(:)).^2)/(Nx*Ny));
hat_x = real(ifft2(ifftshift(Yft_s)));

% Then we correct the MSV output of the L filter in hat_x, keeping the
% previous features with their references values:

Yf = fftshift(fft2(hat_x));
LLLH_desired = orthofeat(5);

t0 = [0 0 0]; %initial guess
fun2 = @compute_f3_n_f4_n_f5_of_normalized_x_BP;
hat_t = fsolve(fun2,t0,options);

% We transform accordingly x
% First, we test that hat_f3 has been normalized:
Yft = Yf.* exp((abs(LH).^2*hat_t(1) + abs(LLH).^2*hat_t(2) + abs(LLLH).^2*hat_t(3)));
Yft_s = Yft/sqrt(mean(abs(Yft(:)).^2)/(Nx*Ny));
hat_x = real(ifft2(ifftshift(Yft_s)));

% Then we correct the MSV output of the L filter in hat_x, keeping the
% previous features with their references values:

Yf = fftshift(fft2(hat_x));
LLH_desired = orthofeat(4);

t0 = [0 0]; %initial guess
fun2 = @compute_f3_n_f4_of_normalized_x_BP;
hat_t = fsolve(fun2,t0,options);

% We transform accordingly x
% First, we test that hat_f3 has been normalized:
Yft = Yf.* exp((abs(LH).^2*hat_t(1) + abs(LLH).^2*hat_t(2)));
Yft_s = Yft/sqrt(mean(abs(Yft(:)).^2)/(Nx*Ny)^2);
hat_x = real(ifft2(ifftshift(Yft_s)));

% Then we correct the MSV output of the L filter in hat_x, keeping the
% previous features with their references values:

Yf = fftshift(fft2(hat_x));
LH_desired = orthofeat(3);

t0 = 0; %initial guess
fun = @compute_f3_of_normalized_x_BP;
hat_t = fsolve(fun,t0,options);

% We transform accordingly x
% First, we test that hat_f3 has been normalized:
Yft = Yf.* exp(abs(LH).^2*hat_t);
Yft_s = Yft/sqrt(mean(abs(Yft(:)).^2)/(Nx*Ny));
hat_x = real(ifft2(ifftshift(Yft_s)));

% Finally, we correct the standard deviation and the mean

w = hat_x*sqrt(orthofeat(2)) + orthofeat(1);


end

