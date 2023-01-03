function f3t = compute_f3_of_normalized_x_BP(t)   
% We apply:
% X(f_i,t) = X(f_i,0) .* exp(2/N |H(f_i)|.^2 t)
% and
% ^t_h(x) = arg_t {f_h(DFT^{-1}{X(f,t)}) = V_h^{ref}}.
% In addition, we impose the Euclidean norm normalization, so the
% adjustment is double: in the pixel domain and at the output of the
% Low-pass filter L.

global Yf
global LH_desired
global LH

N = size(Yf,1);
Yft = Yf.* exp(abs(LH).^2*t);
Yft_s = Yft/sqrt(mean(abs(Yft(:)).^2)/N^2);
f3t = mean(abs(LH(:).*Yft_s(:)).^2/N^2) - LH_desired;
end



