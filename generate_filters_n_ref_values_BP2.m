function [LH,LLH,LLLH,LLLL,LH_ref,LLH_ref,LLLH_ref,LLLL_ref] = generate_filters_n_ref_values_BP2(SizeIm)
% This function computes the spectral responses of a set of three band-pass filters
% (LH, LLH, LLLH) and of a very low-pass filter (LLLL).
% It also computes the mean square expected values at the output of these
% filters (LH_ref, etc.), when applying at the input a
% white standardized noise input of size Ny x Nx.
%
% INPUT:
%   SizeIm = size(im) = [Ny Nx]:  image size (Ny x Nx)
%
% OUTPUT:
%   LH: Centered (with fftshift) response of the first band-pass filter
%   LLH: Centered (with fftshift) response of the (very low) band-pass filter
%   LLLH: Centered (with fftshift) response of the (ultra low) band-pass filter
%   LLLL: Centered (with fftshift) response of the (extreme low)-pass filter
%   LH_ref: MSV reference value at the output of LH
%   LLH_ref: MSV reference value at the output of LLH
%   LLLH_ref: MSV reference value at the output of LLLH
%   LLLL_ref: MSV reference value at the output of LLLL
%
% These filters, together with their complementary in square value, are a
% Parseval frame, which means that MSV_pixel = MSV_H + MSV_L, and MSV_L =
% MSV_LH + MSV_LL, etc., where the LH filter is a band pass filter:
% L.*H1 = L.*sqrt(1-L1.^2), and LL = L.*L1, being:
% L = cos(pi*f).*(f<=1/2), and L1 = cos(2*pi*f).*(f<=1/4), and f the radial
% spatial frequency.

% JPM.

Ny = SizeIm(1); Nx = SizeIm(2); % Watch out: they are assumed to be even, for now
[u,v] = meshgrid(-Nx/2:Nx/2-1,-Ny/2:Ny/2-1);
u = u/Nx; v = v/Ny;

f = sqrt(u.^2 + v.^2);

% First scale filter:
L0 = cos(pi*f).*(f<=1/2);
H0 = sqrt(1 - L0.^2); % High-pass filter (implicit)
L = L0; % First low pass filter  (implicit)
% Second scale and beyond are band-limited.
L1 = cos(2*pi*f).*(f<=1/4); 
H1 = sqrt(1 - L1.^2);
LH = L.*H1;  % Band-pass filter
LL = L.*L1; % Second low pass filter (implicit)
% Third scale
L2 = cos(4*pi*f).*(f<=1/8); 
H2 = sqrt(1 - L2.^2);
LLH = LL.*H2;  % Band-pass filter
LLL = LL.*L2; % Third low pass filter (implicit)
% Fourth scale
L3 = cos(8*pi*f).*(f<=1/16); 
H3 = sqrt(1 - L3.^2);
LLLH = LLL.*H2;  % Band-pass filter
LLLL = LLL.*L2; % Fourth low pass filter

% v_2^{ref} = sigma^2 = 1; v_1^{ref} = mu = 0
LH_ref = sum(L(:).^2)/(Ny*Nx); % = v_3^{ref};  
LLH_ref = sum(LL(:).^2)/(Ny*Nx); % = v_4^{ref};
LLLH_ref = sum(LLL(:).^2)/(Ny*Nx); % = v_4^{ref};
LLLL_ref = sum(LLLL(:).^2)/(Ny*Nx); % = v_4^{ref};

end
