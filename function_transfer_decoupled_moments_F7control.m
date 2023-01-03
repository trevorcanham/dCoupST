% Compute kernel for diffusion transfer
% JPM
function [kernel, source_features, reference_features] = function_transfer_decoupled_moments_F7control(S, R, Q)

    global Xfrec Kern

    global LHs LLHs LLLHs LLLLs
    global LHr LLHr LLLHr LLLLr

    reference_features = [];
    source_features = [];

    [Sh,Sw,Sc] = size(S);
    [Rh,Rw,Rc] = size(R);    
    
%% Set parameters for preprocessing
% Quantization

    Sf = reshape(S,size(S,1)*size(S,2),size(S,3));
    Rf = reshape(R,size(R,1)*size(R,2),size(R,3));

    Sf = perturbate_for_moments(Sf,Q);
  
    % Undo gamma compression (to obtain values proportional to physical radiance)
    mean_S = mean(Sf).';
    mean_R = mean(Rf).';

    % The first set of features are the mean R, G and B radiance values
    reference_features = [reference_features; mean_R];
    source_features = [source_features; mean_S];

    S2f = Sf;
    S2 = reshape(S2f,size(S,1),size(S,2),size(S,3));
    R2 = reshape(Rf,size(R,1),size(R,2),size(R,3));

%% Compute kernel
    % Now it imposes the spectral features from R to S
    % It does it by computing the PSF that fixes the decoupled variance at a set
    % of band-pass filters. Then, it applies that PSF to each of the R, G and B
    % channesl.

    L_S20 = S2;
    L_R0 = R2;

    im_old = L_S20;

    % Here I measure (in the reference L_R) and impose (to the source L_S2) the
    % orthovariances at the output of 3 band-pass isotropic filters. Because
    % these measurements are subject to zero mean and
    % unit variance, they actually control also the variance at the output of 
    % implicit high frequency and low frequency filters.

    % Now we extract the orthofeatures of the reference image
    [orthofilterR] = normalize_n_measure_orthofilters_BP2(L_R0,LHr,LLHr,LLLHr,LLLLr);
    % And now we use the same function, but for normalizing the source
    [orthofilterS,hat_x] = normalize_n_measure_orthofilters_BP2(L_S20,LHs,LLHs,LLLHs,LLLLs);

    % Finally, we impose the orthofeatures measured in the origin to the source
    [L_S20] = impose_orthofilters_BP2(hat_x,orthofilterR,LHs,LLHs,LLLHs,LLLLs);

    % Now we compute the equivalent filter adjusting the spectral properties
    im = L_S20;
    Im_old = fft2(im_old);
    Im = fft2(im);
    % Here it is important to realize that the transformed image is a filtered version of the original 
    lambda = 1e-8;
    Heq = Im./(Im_old + lambda);
    heq = fftshift(real(ifft2(Heq)));

    % We need to crop here the central part
    L = 25;
    kernel = heq((Sh/2+1)-L:(Sh/2+1)+L,(Sw/2+1)-L:(Sw/2+1)+L);

end

