function [cart, mask] = MRI_Cartesian(img, nlines, npoints, mType, mPercent) 
    len = length(img);
    n = [len/nlines, len/npoints];
    m = floor(len*n);
    img_zeros = zeros(m(1), m(2));
    img_zeros(1:len, 1:len) = img;
    kspace = fftshift(fft2(img_zeros)); 
    grid = zeros(m(1),m(2));

    % sample intervals
    sample = interp2(kspace, (m(2)/2-len/2:n(2):m(2)/2+len/2-1)',(m(1)/2-len/2:n(1):m(1)/2+len/2-1));
    sample_size = size(sample);
    sample_zeros = zeros(sample_size(1), sample_size(2));
    
    mask = getMask(sample_size, mType, mPercent);
    
    for i=1:1:sample_size(1)
        for j=1:1:sample_size(2)
            if mask(i,j) == 1
                sample_zeros(i, j) = sample(i,j);
            end
        end
    end
    
    grid(m(1)/2-sample_size(1)/2+1:(m(1)/2+sample_size(1)/2),  m(2)/2-sample_size(2)/2+1:(m(2)/2+sample_size(2)/2)) = sample_zeros;
    grid(isnan(grid)) = 0;
    
    % Fourier transform
    X = (ifft2(fftshift(grid)));
    X = abs((X));
    
    res_X = imresize(X, size(X)./n); 
    % rescale img
    cart = res_X;
    
    cart = cart/(max(cart(:))) * 255;

    
end
