% apply weighted average between three color channels
% Copyright Trevor D. Canham 2022
function out = treeChanWeight(im1,im2,weights)

    out = im1;
    
    for i = 1:size(im1,3)   
        out(:,:,i) = im1(:,:,i) .* weights(i) + im2(:,:,i) .* (1-weights(i));
    end
    
end