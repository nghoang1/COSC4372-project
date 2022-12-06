function [outputImg, ourMask] = MRI_Radial(img, lines, pointsOnLine, maskType, maskPercent)
    imageSize = size(img);
    sample = (imageSize(1)/pointsOnLine);
    Nsize = sample * imageSize(1) * 3;

    imageGrid = zeros(Nsize, Nsize);
    imageGrid(1:imageSize(1), 1:imageSize(1)) = img;
    kSpace = fftshift(fft2(imageGrid));

    i=1;
    j=1;
    for kSpace =- Nsize/2:sample:Nsize/2
       for theta = 0 : pi/lines : (pi-pi / lines)

           radial_X(i, j) = kSpace * cos(-theta) + Nsize / 2;
           radial_Y(i, j) = kSpace * sin(-theta) + Nsize / 2;
           
           i = i + 1;
       end
       j = j + 1;
       i = 1;
    end
    
    %radial view sampling
    radialView = interp2(F, radial_X, radial_Y, 'bicubic');
    radialView(isnan(radialView)) = 0;

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

    %Inverse Fourier Transformation
    IR = zeros(size(maskedRadialView));    
    IR = InverseFourierTransformation(IR, maskedRadialView, lines);

    %image recreation
    imgCreator = iradon(IR', 180/delT);
    imgCreator = rot90(imgCreator(16:sizeOfImage+15, 16:sizeOfImage+15),2);
    outputImg = imgCreator;  
end

function [radialView] = InverseFourierTransformation(radialView, maskedView, lines)
    for i = 1:lines
        radialView = fliplr(fftshift((abs(ifft((maskedView(i, :)))))));
    end
end