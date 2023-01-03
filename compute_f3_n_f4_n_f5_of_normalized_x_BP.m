function f3t_n_f4t_n_f5t = compute_f3_n_f4_n_f5_of_normalized_x_BP(t)   
% We apply:
% X(f_i,t1,t2,t3) = X(f_i,0,0,0) .* exp(2/N*( |H1(f_i)|.^2 t1 + |H2(f_i)|.^2 t2 + |H3(f_i)|.^2 t3))
% and
% ^t_h(x) = arg_(t1,t2,t3) {f_h(DFT^{-1}{X(f,t1,t2,t3)}) = V_h^{ref} = [V_1^{ref},V_2^{ref},V_3^{ref}].
% In addition, we impose the Euclidean norm normalization, so the
% adjustment is quadruple: in the pixel domain and at the output of the
% three filters, H1, H2 and H3.

global Yf
global LH_ref
global LLH_ref
global LLLH_desired
global LH
global LLH
global LLLH

N = size(Yf,1);
Yft = Yf.* exp((abs(LH).^2*t(1) + abs(LLH).^2*t(2) + abs(LLLH).^2*t(3)));
Yft_s = Yft/sqrt(mean(abs(Yft(:)).^2)/N^2);
f3t = mean(abs(LH(:).*Yft_s(:)).^2/N^2) - LH_ref;
f4t = mean(abs(LLH(:).*Yft_s(:)).^2/N^2) - LLH_ref;
f5t = mean(abs(LLLH(:).*Yft_s(:)).^2/N^2) - LLLH_desired;
f3t_n_f4t_n_f5t = [f3t f4t f5t];
end