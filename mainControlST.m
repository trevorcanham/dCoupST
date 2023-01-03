% Using Decoupled Features for Photorealistic Style Transfer
% Arxiv Preprint: https://doi.org/10.48550/arXiv.2212.02953
% Given source and target images (S and T), a transform to match decoupled visual response features is derived in the
% form of a 3DLUT and a convolution kernel and applied to a working image
% (W).
% Check readme for inputs/outputs, or run dCoupSTplus to set parameters, call images
% Copyright Trevor D. Canham, Adrián Martín, Marcelo Bertalmío, and Javier Portilla 2022

function [out,moments,kernel,lut] = mainControlST(W,S,R,params)

    addpath('utility');

    global LHs LLHs LLLHs LLLLs % init. filters for diffusion transfer
    global LHr LLHr LLLHr LLLLr
    global LH_ref LLH_ref LLLH_ref
    
    maxVal = (2^params.depth)-1; % init. input image max code value
    
    % read images
    W = double(W)./maxVal; 
    W = imresize(W,params.resizeF);
    W(W<0) = 0;
    S = double(S)./maxVal;  
    S = imresize(S,params.resizeF);
    S(S<0) = 0;    
    R = double(R)./maxVal;        
    R = imresize(R,params.resizeF);
    R(R<0) = 0;  
    [Sh,Sw,Sr] = size(S);
    [Wh,Ww,Wr] = size(W);
    [Rh0,Rw0,Rr0] = size(R);
    
    % Conform all inputs to n X m X 3 data structure
    if Sr < 3
        S = cat(3,S(:,:,1),S(:,:,1),S(:,:,1));
        Sr = 3;
    end
    
    if Wr < 3
        W = cat(3,W(:,:,1),W(:,:,1),W(:,:,1));
        Wr = 3;
    end

    if Rr0 < 3
        R = cat(3,R(:,:,1),R(:,:,1),R(:,:,1));
        Rr0 = 3;
    end
          
    LoadGamuts; % import primary matrices, assign parameters
    
    switch params.sCSin{:}
        case 'sRGB'
            params.sCSmat = M_srgb;
        case 'BT.2020'
            params.sCSmat = M_BT2020;
        case 'ProPhoto'
            params.sCSmat = M_ProPhoto;
        case 'ARRI'
            params.sCSmat = M_ARRI;
    end
   
    switch params.rCSin{:}
        case 'sRGB'
            params.rCSmat = M_srgb;
        case 'BT.2020'
            params.rCSmat = M_BT2020;
        case 'ProPhoto'
            params.rCSmat = M_ProPhoto;
        case 'ARRI'
            params.rCSmat = M_ARRI;
    end    
    
    params.csOutMat = M_srgb;
    switch params.csOut{:}
        case 'sRGB'
            params.csOutMat = M_srgb;
        case 'BT.2020'
            params.csOutMat = M_BT2020;
        case 'DCI.P3.D65'
            params.csOutMat = M_dcip3_D65;
        case 'DCI.P3.D50'
            params.csOutMat = M_dcip3;
    end     
    
    
    % if not already computed .. compute LUT and kernel
    if ~isa(params.moments,'cell')
        
       % force source and reference into equal, square dimensions for
       % kernel computation
        Rre = imresize(R,[Sh,Sw]);
        [Rh,Rw,Rr] = size(Rre);
        N = min([Sh Sw Rh Rw]);
        N = 2*floor(N/2);
        
        if Sh/2~=round(Sh/2)
            S = S(1:end-1,:,:);
            Sh = Sh -1;
        end
        
        if Sw/2~=round(Sw/2)
            S = S(:,1:end-1,:);
            Sw = Sw - 1;
        end    
        
        if Rh/2~=round(Rh/2)
            Rre = Rre(1:end-1,:,:);
            Rh = Rh - 1;
        end    
        
        if Rw/2~=round(Rw/2)
            Rre = Rre(:,1:end-1,:);
            Rw = Rw - 1;
        end  
        
        Scrop = S;
        Rcrop = Rre;
        
        % convert to RGB display linear 
        ScropLin = ((Scrop).^params.sGamIn);
        RcropLin = ((Rcrop).^params.rGamIn);        
        
        sKernel = ScropLin;
        
        % pull luminance for kernel computation
        ScropLinChan =  0.2126.*ScropLin(:,:,1) + 0.7152.*ScropLin(:,:,2) + 0.0722.*ScropLin(:,:,3);
        RcropLinChan =  0.2126.*RcropLin(:,:,1) + 0.7152.*RcropLin(:,:,2) + 0.0722.*RcropLin(:,:,3);
        
        % Create the corresponding filters and reference values
        if min(Sw,Sh)>min(Rw,Rh) % it chooses the largest image for computing the common reference values
            [LHs,LLHs,LLLHs, LLLLs,LH_ref,LLH_ref,LLLH_ref,LLLL_ref] = generate_filters_n_ref_values_BP2([Sh Sw]);
            [LHr,LLHr,LLLHr, LLLLr] = generate_filters_n_ref_values_BP2([Rh Rw]);
        else
            [LHr,LLHr,LLLHr, LLLLr,LH_ref,LLH_ref,LLLH_ref,LLLL_ref] = generate_filters_n_ref_values_BP2([Rh Rw]);    
            [LHs,LLHs,LLLHs, LLLLs] = generate_filters_n_ref_values_BP2([Sh Sw]);
        end

        % compute kernel based on source and reference 
        [params.kernel]  = function_transfer_decoupled_moments_F7control(ScropLinChan, RcropLinChan, maxVal);       
        params.moments = cell(Sr,2);
        
        % initialize LUT
        [a,b,c] = meshgrid(linspace(0,1,params.lutSize));
        idLut = horzcat(b(:),a(:),c(:));                
        % convert to display linear RGB
        sLin = (S.^params.sGamIn);
        rLin = (R.^params.sGamIn);  
        lutLin = (idLut.^params.sGamIn);
        sLinColor = sLin;
        lutLinColor = lutLin;
        
        % transfer mean on RGB channels, if color channels are enabled
        if sum(params.moms(1,:)) == 3       
            docolor = 1;
            for i = 1:3
                sChan = sLin(:,:,i);
                rChan = rLin(:,:,i);
                sLinColor(:,:,i) = sLin(:,:,i).*(mean(rChan(:))/mean(sChan(:)));
                lutLinColor(:,i) = lutLin(:,i).*(mean(rChan(:))/mean(sChan(:)));
            end
            sLinColor = min(max(sLinColor,0),1);
            lutLinColor = min(max(lutLinColor,0),1);
        else
            docolor = 0;
            sLinColor = min(max(sLinColor,0),1);
            lutLinColor = min(max(lutLinColor,0),1);
        end   

        ST1 = sLinColor;
        
        % apply kernel
        for i = 1:Sr
            
            if params.moms(5,i) > 0
                ST1(:,:,i) = conv2(sLinColor(:,:,i),params.kernel,'same');
            else
                ST1(:,:,i) = sLinColor(:,:,i);
             end
                
        end

        % weight vs. original (diffusion weights can be set anywhere in 0-1 range for a partial transfer)
        sLinColor = treeChanWeight(ST1,sLinColor,params.moms(5,:));        
        
        % convert to IPT for color transfer
        if strcmp(params.csInt{:},'IPT')

            sLinColorVect = reshape(sLinColor,Sh*Sw,Sr);
            rLinVect = reshape(rLin,Rh0*Rw0,Rr0);
            sXYZ = ConvertRGBtoXYZ(sLinColorVect,0);
            lutXYZ = ConvertRGBtoXYZ(lutLinColor,0);
            rXYZ = ConvertRGBtoXYZ(rLinVect,0);
            sIntVect = ConvertXYZtoIPT(sXYZ,0);
            rIntVect = ConvertXYZtoIPT(rXYZ,0);
            lutInt = ConvertXYZtoIPT(lutXYZ,0);
            ST1weight = reshape(sIntVect,Sh,Sw,Sr);
            rInt = reshape(rIntVect,Rh0,Rw0,Rr0);        
            
        else % RGB case
            
            rInt = rLin;
            ST1weight = sLinColor;
            lutInt = lutLinColor;
        
        end
        
        % apply moments according to settings
        ST2 = ST1weight;
        lutST = lutInt;
        
        for i = 1:Sr
        
            if params.moms(1,i) < 1
                ST2(:,:,i) = ST1weight(:,:,i);
                source_features = zeros(4,1);
                reference_features = zeros(4,1);
                
            elseif params.moms(2,i) < 1       
                [ST2(:,:,i), source_features, reference_features,lutST(:,i)] = toggle(ST1weight(:,:,i), rInt(:,:,i), 1, maxVal, lutInt(:,i));
            elseif params.moms(3,i) < 1       
                [ST2(:,:,i), source_features, reference_features,lutST(:,i)] = toggle(ST1weight(:,:,i), rInt(:,:,i), 2, maxVal, lutInt(:,i));
            elseif params.moms(4,i) < 1       
                [ST2(:,:,i), source_features, reference_features,lutST(:,i)] = toggle(ST1weight(:,:,i), rInt(:,:,i), 3, maxVal, lutInt(:,i));
            else
                [ST2(:,:,i), source_features, reference_features,lutST(:,i)] = toggle(ST1weight(:,:,i), rInt(:,:,i), 4, maxVal, lutInt(:,i));
            end
        
            params.moments{i,1} = source_features;
            params.moments{i,2} = reference_features;
            
        end
        
        % convert back from IPT
        if strcmp(params.csInt{:},'IPT')
            
            outLinVect = reshape(ST2,Sh*Sw,Sr);
            outXYZ = ConvertXYZtoIPT(outLinVect,1);
            lutXYZ = ConvertXYZtoIPT(lutST,1);       
            outLinVect = ConvertRGBtoXYZ(outXYZ,1);
            lutOut = ConvertRGBtoXYZ(lutXYZ,1);
            outLin = reshape(outLinVect,Sh,Sw,Sr);  
            
        else 
            
            outLin = ST2;
            lutOut = lutST;
        
        end   

        % limit to 0-1 range
        outLin = min(max(outLin,0),1);
        lutOut = min(max(lutOut,0),1);
        outLinColor = outLin;
        lutOutColor = lutOut;
        
        % RGB means post-adjustment
        if docolor 
            for i = 1:3
                outChan = outLin(:,:,i);
                rChan = rLin(:,:,i);
                outLinColor(:,:,i) = outLin(:,:,i).*(mean(rChan(:))/mean(outChan(:)));
                lutOutColor(:,i) = lutOut(:,i).*(mean(rChan(:))/mean(outChan(:)));
            end
            outLinColor = min(max(outLinColor,0),1);
            %lutOutColor = min(max(lutOutColor,0),1);
        else  
            outLinColor = min(max(outLinColor,0),1);
            %lutOutColor = min(max(lutOutColor,0),1);
        end    
        
        % convert to display-ready representation
        outS = real((outLinColor).^(1/params.gamOut)); % .. for curiosity
        lut = real((lutOutColor).^(1/params.gamOut)); 
        %lut = min(max(lut,0),1);

        params.luts = lut;
                       
    end

    
    % apply LUT, kernel to working image for output
    % convert S,R to display linear RGB   
    sLin = (W.^params.sGamIn); 
    ST1 = sLin;

    % apply auto covariance
    for i = 1:Sr

        if params.moms(5,i) > 0
            ST1(:,:,i) = conv2(sLin(:,:,i),params.kernel,'same');
        else
            ST1(:,:,i) = sLin(:,:,i);
        end

    end

    % weight vs. original
    sLinKern = treeChanWeight(ST1,sLin,params.moms(5,:));        

    % return to input encoding,  limit to real numbers, 0-1 range
    sLinKern = min(max(sLinKern,0),1);
    sKern = real((sLinKern).^(1/params.gamOut));

    % apply color transform 3DLUT
    out = imlut(sKern,params.luts,'3D','standard','BGR');
    
    moments = params.moments;
    kernel = params.kernel;   
    lut = params.luts;    
    
end