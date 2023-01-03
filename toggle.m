% compute and apply nested moments for color transfer to source image and LUT
% Adapted from JPM code to control which moments to transfer
% Since the transfer is nested, users can select to transfer 1st moment,
% 1st and 2nd moment, 1st,2nd, and 3rd moments or 1st,2nd,3rd, and
% 4th moments. The user may not, however, transfer only higher moments
% without transferring all of the lower ones.
% For 1st, 2nd and 3rd moments the color transfer is analytical and can be
% easily applied to the lookup table, however the 4th moment color transfer
% is approximated with a spline function, which is later applied to the
% LUT.
% Copyright Trevor D. Canham 2022
function [out, source_features, reference_features, lut] = toggle(x1, x2, Nmoms, maxVal, idLut)
    
    source_features = zeros(4,1);
    reference_features = zeros(4,1);
    
    Q = maxVal+1;

    x1pert = perturbate_for_moments(x1(:),Q);
    x1 = reshape(real(x1pert),size(x1,1),size(x1,2));

    switch Nmoms
        
        case 1
            
            % transfer reference mean to source
            mu1 = mean(x1(:));
            mu2 = mean(x2(:)); 
            xNorm1 = x1 - mu1;
            lutNorm1 = idLut - mu1; 
            out = xNorm1 + mu2;
            lut = lutNorm1 + mu2;
            source_features(1,:) = mu1;
            reference_features(1,:) = mu2;
            
        case 2
            
            % transfer reference mean and variance to source
            mu1 = mean(x1(:));
            mu2 = mean(x2(:));   
            xNorm1 = x1 - mu1;
            lutNorm1 = idLut - mu1;
            var1 = mean(xNorm1(:).^2);
            xNorm2 = x2 - mu2;
            var2 = mean(xNorm2(:).^2);
            out = (xNorm1.*sqrt(var2/var1))+mu2; 
            lut = (lutNorm1.*sqrt(var2/var1))+mu2; 
            source_features(1:2,:) = [mu1,var1];
            reference_features(1:2,:) = [mu2,var2];
             
        case 3
            
            % transfer reference mean, variance, and skew to source
            mu1 = mean(x1(:));
            xNorm1 = x1 - mu1;
            lutNorm1 = idLut - mu1;
            var1 = mean(xNorm1(:).^2);
            xVarNorm1 = xNorm1./sqrt(var1);
            lutVarNorm1 = lutNorm1./sqrt(var1);
            sk1 = mean(xVarNorm1(:).^3);
            mu2 = mean(x2(:)); 
            xNorm2 = x2 - mu2;
            var2 = mean(xNorm2(:).^2);
            xVarNorm2 = xNorm2./sqrt(var2);
            sk2 = mean(xVarNorm2(:).^3);
            %impose moments
            if sk1 > sk2
                t_guess = [(1-1e-8)/min(xVarNorm1(:)) 0];
            else    
                t_guess = [0 (1-1e-8)/max(xVarNorm1(:))];
            end
            t_sk2 = fzero(@(t) skew_adj(t,real(xVarNorm1(:)))-sk2,real(t_guess));
            % Impose skew
            xSkew = xVarNorm1./(1-t_sk2*xVarNorm1);
            lutSkew = lutVarNorm1./(1-t_sk2*lutVarNorm1);
            % Impose mean and variance
            xSkewMuNorm = xSkew - mean(xSkew(:));
            lutSkewMuNorm = lutSkew - mean(xSkew(:));
            xSkewMuNormVar = mean(xSkewMuNorm(:).^2);
            out = xSkewMuNorm.*sqrt(var2/xSkewMuNormVar)+mu2;
            lut = lutSkewMuNorm.*sqrt(var2/xSkewMuNormVar)+mu2;
            source_features(1:3,:) = [mu1,var1,sk1];
            reference_features(1:3,:) = [mu2,var2,sk2];            

        case 4
            
            % transfer reference mean, variance, skew, and ortho-kurtosis
            % to source
            mu1 = mean(x1(:));
            x1Norm1 = x1 - mu1;
            lut1Norm1 = idLut - mu1;
            var1 = mean(x1Norm1(:).^2);
            x1Norm2 = x1Norm1./sqrt(var1);
            lut1Norm2 = lut1Norm1./sqrt(var1);
            sk1 = mean(x1Norm2(:).^3);
            [ok1, x1Norm3, ~, lut1Norm3] = fastOrthoKurtosisLut(real(x1Norm2),lut1Norm2);
            source_features = [mu1,var1,sk1,ok1];

            mu2 = mean(x2(:));
            x2Norm1 = x2 - mu2;
            var2 = mean(x2Norm1(:).^2);
            x2Norm2 = x2Norm1./sqrt(var2);
            sk2 = mean(x2Norm2(:).^3);
            [ok2, ~] = fastOrthoKurtosis(real(x2Norm2));
            reference_features = [mu2,var2,sk2,ok2];

            % Impose the decoupled moments to x1Norm3 in reverse order (from high
            % order to low order)
            rel_tol = 1e-6;
            [~,~,out] = impose_n_measure_features_reverse(x1Norm3,reference_features,rel_tol);
            lut = spline(x1Norm3(:),out(:),lut1Norm3);
            
    end

    out = real(out);
    lut = real(lut');

end



