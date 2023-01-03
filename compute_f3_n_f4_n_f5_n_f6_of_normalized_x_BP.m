function f3t_n_f4t_n_f5t_n_f6t = compute_f3_n_f4_n_f5_n_f6_of_normalized_x_BP(t)   
% We apply:
% X(f_i,t1,t2,t3,t4) = X(f_i,0,0,0,0) .* exp(2/N*( |H1(f_i)|.^2 t1 + |H2(f_i)|.^2 t2 + |H3(f_i)|.^2 t3 + |H4(f_i)|.^2 t4))
% and
% ^t_h(x) = arg_(t1,t2,t3,t4) {f_h(DFT^{-1}{X(f,t1,t2,t3,t4)}) = V_h^{ref} = [V_1^{ref},V_2^{ref},V_3^{ref},V_4^{ref}].
% In addition, we impose the Euclidean norm normalization, so the
% adjustment is 5 fold: in the pixel domain and at the output of the
% three filters, H1, H2, H3 and H4.

global Yf
global LH_ref
global LLH_ref
global LLLH_ref
global LLLL_desired
global LH
global LLH
global LLLH
global LLLL

N = size(Yf,1);
Yft = Yf.* exp((abs(LH).^2*t(1) + abs(LLH).^2*t(2) + abs(LLLH).^2*t(3) + abs(LLLL).^2*t(4)));
Yft_s = Yft/sqrt(mean(abs(Yft(:)).^2)/N^2);
f3t = mean(abs(LH(:).*Yft_s(:)).^2/N^2) - LH_ref;
f4t = mean(abs(LLH(:).*Yft_s(:)).^2/N^2) - LLH_ref;
f5t = mean(abs(LLLH(:).*Yft_s(:)).^2/N^2) - LLLH_ref;
f6t = mean(abs(LLLL(:).*Yft_s(:)).^2/N^2) - LLLL_desired;
f3t_n_f4t_n_f5t_n_f6t = [f3t f4t f5t f6t];
end

