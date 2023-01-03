function [orthofeat,hat_x,feat] = normalize_n_measure_orthofilters_BP2(x,LHr,LLHr,LLLHr,LLLLr) %,L,LL,L_ref)
% This function computes the orthofeatures from the following hierarchical
% set:
% 1) the mean
% 2) MSV (mean Square value) in the pixel domain
% 3) MSV at the output of a high-frequency band-pass filter, LH
% 4) MSV at the output of a low-frequency band-pass filter LLH
% 5) MSV at the output of a very low-frequency band-pass filter, LLLH
% 6) MSV at the output of a very low-pass filter LLLL
%
% INPUT:
%   x:  image (gray level)
% (Implicitely, through global variables)
%   LH:  spectral (centered) response of a band-pass filter
%   LLH:  spectral (centered) response of a low-frequency band-pass filter
%   LLLH:  spectral (centered) response of a very low-frequency band-pass filter
%   LLLL:  spectral (centered) response of a very low-pass filter
%   LH_ref: reference MSV at the output of LH
%   LLH_ref = reference MSV at the output of LLH
%   LLLH_ref = reference MSV at the output of LLLH
%
% OUTPUT:
%   orthofeat: corresponding set of orthofeatures
%   hat_x: normalized x, up to the N-1 (fifth) feature.
%   feat: features
%
% (the last feature does not need to be normalized, as it is going to be
% changed later on)

global Yf
global LH_ref
global LLH_ref
global LLLH_ref
global LH_desired
global LLH_desired
global LLLH_desired

global LH
global LLH
global LLLH
global LLLL

LH = LHr;
LLH = LLHr;
LLLH = LLLHr;
LLLL = LLLLr;

[Ny,Nx] = size(x);
if (2*round(Ny/2)~=Ny)||(2*round(Nx/2)~=Nx)
    error('Input image must have even size dimensions');
    return;
end

% First the mean and variance (standardize)
f_1 = mean(x(:));
hat_f_1 = f_1;
hat_x_1 = x - f_1;
f_2 = mean(x(:).^2);
hat_f_2 = mean(hat_x_1(:).^2); % ^f_2(x) = f_2(^x_1) 
hat_x_2 = hat_x_1/sqrt(hat_f_2);

% And measure the third feature and orthofeature
Xf = fftshift(fft2(x));
xLH = real(ifftshift(ifft2(Xf.*LH)));
f_3 = mean(xLH(:).^2);
y = hat_x_2;
Yf = fftshift(fft2(y));
yLH = real(ifftshift(ifft2(Yf.*LH)));
hat_f_3 = mean(yLH(:).^2);  % ^f_3(x) = f_3(^x_2) 

% Now we normalize also the MSV at the output of the L filter (low
% frequency), and obtain the fourth feature and orthofeature
LH_desired = LH_ref;
t0 = 0; %initial guess
fun = @compute_f3_of_normalized_x_BP;
%options = optimoptions('fsolve','FunctionTolerance',1.0e-16,'OptimalityTolerance',1e-9);
options = optimoptions('fsolve','FunctionTolerance',1.0e-12,'OptimalityTolerance',1e-6);
hat_t = fsolve(fun,t0,options);

% First, we test that hat_f3 has been normalized:
% error = compute_f3_of_normalized_x(hat_t);
Yft = Yf.* exp(abs(LH).^2*hat_t);
Yft_s = Yft/sqrt(mean(abs(Yft(:)).^2)/(Ny*Nx));
hat_x_3 = real(ifft2(ifftshift(Yft_s)));

% Now it measures the fourth feature and orthofeature
xLLH = real(ifftshift(ifft2(Xf.*LLH)));
f_4 = mean(xLLH(:).^2);
z = hat_x_3;
Zf = fftshift(fft2(z));
zLLH = real(ifft2(ifftshift(Zf.*LLH)));
hat_f_4 = mean(zLLH(:).^2);  % ^f_4(x) = f_4(^x_3) 

% Now we normalize also the MSV at the output of the LL filter (very low
% frequency), and obtain the fifth feature and orthofeature
t0 = [0 0]; %initial guess
LLH_desired = LLH_ref;
fun2 = @compute_f3_n_f4_of_normalized_x_BP;
hat_t = fsolve(fun2,t0,options);

% We transform accordingly x
%
%error = compute_f3_n_f4_of_normalized_x(hat_t)
Yft = Yf.* exp((abs(LH).^2*hat_t(1) + abs(LLH).^2*hat_t(2)));
Yft_s = Yft/sqrt(mean(abs(Yft(:)).^2)/(Ny*Nx));
hat_x_4 = real(ifft2(ifftshift(Yft_s)));

% Now it measures the fifth feature and orthofeature
xLLLH = real(ifftshift(ifft2(Xf.*LLLH)));
f_5 = mean(xLLLH(:).^2);
z = hat_x_4;
Zf = fftshift(fft2(z));
zLLLH = real(ifft2(ifftshift(Zf.*LLLH)));
hat_f_5 = mean(zLLLH(:).^2);  % ^f_5(x) = f_5(^x_4) 

% Now we normalize also the MSV at the output of the LL filter (very low
% frequency), and obtain the fifth feature and orthofeature
t0 = [0 0 0]; %initial guess
LLLH_desired = LLLH_ref;
fun2 = @compute_f3_n_f4_n_f5_of_normalized_x_BP;
hat_t = fsolve(fun2,t0,options);

% We transform accordingly x
%
%error = compute_f3_n_f4_of_normalized_x(hat_t)
Yft = Yf.* exp((abs(LH).^2*hat_t(1) + abs(LLH).^2*hat_t(2) + abs(LLLH).^2*hat_t(3)));
Yft_s = Yft/sqrt(mean(abs(Yft(:)).^2)/(Ny*Nx));
hat_x_5 = real(ifft2(ifftshift(Yft_s)));

% Now it measures the sixth feature and orthofeature
xLLLL = real(ifftshift(ifft2(Xf.*LLLL)));
f_6 = mean(xLLLL(:).^2);
z = hat_x_5;
Zf = fftshift(fft2(z));
zLLLL = real(ifft2(ifftshift(Zf.*LLLL)));
hat_f_6 = mean(zLLLL(:).^2);  % ^f_6(x) = f_6(^x_5) 


feat = [f_1 f_2 f_3 f_4 f_5 f_6];
orthofeat = [hat_f_1 hat_f_2 hat_f_3 hat_f_4 hat_f_5 hat_f_6];

hat_x = hat_x_5; % No need to compute hat_x_6


end

