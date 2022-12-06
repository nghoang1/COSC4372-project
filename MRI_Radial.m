function [finalImg, ourMask] = MRI_Radial(img, lines, ppl, maskType, maskPercent)
  
    
    sizeOfImage = size(img);
    
    
    sample = (sizeOfImage(1)/ppl);
    
    Nsize = sizeOfImage(1) * 3 * sample;

    ourArr = zeros(Nsize, Nsize);
    ourArr(1:sizeOfImage(1), 1:sizeOfImage(1)) = img;
    kSpaceConvert = fftshift(fft2(ourArr));
    
    
    i=1;
    j=1;
    delT = lines;

    for kSpaceConvert=-Nsize/2:sample:Nsize/2

       for theta = 0:pi/delT:(pi-pi/delT)

           radialX(i, j) = kSpaceConvert*cos(-theta)+ Nsize/2;
           radialY(i, j) = kSpaceConvert*sin(-theta)+ Nsize/2;
           
           i = i+1;
       end

       j = j+1;
       i = 1;

    end
    
    
    % Start the sampling for the radial view

    radialView = interp2(fftshift(fft2(ourArr)), radialX, radialY, 'bicubic');
    radialView(isnan(radialView)) = 0;
    %interpolate

    S = size(radialView);
    maskedRadialView = zeros(S);
    ourMask = getMask(S, maskType, maskPercent);
    for i=1:1:S(1)
        for j=1:1:S(2)
            if ourMask(i,j) == 1
                maskedRadialView(i, j) = radialView(i,j);
            end
        end
    end

    
    IR = zeros(size(maskedRadialView));
    
    %Inverse fast fourier transform

    for i = 1:delT

       IR(i, :) =fliplr(fftshift((abs(ifft((maskedRadialView(i, :))))))); 

    end
    %recreating the image for display
    imgCreator = iradon(IR', 180/delT);
    imgCreator = rot90(imgCreator(16:sizeOfImage+15, 16:sizeOfImage+15),2);

    finalImg = imgCreator;
       
end
