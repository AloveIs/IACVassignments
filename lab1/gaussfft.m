function out = gaussfft(img,t)

    %[X,Y] = meshgrid(1:size(img,1),1:size(img,1));

    [X,Y] = meshgrid(-size(img,1)/2:(size(img,1)/2-1),-size(img,1)/2:(size(img,1)/2-1));
    Z = mvnpdf([X(:) Y(:)],0,eye(2,2)*t);
    
    Z = reshape(Z,size(X,1),size(X,1));
    %figure();
    %surf(X,Y,Z);
    Zhat = fft2(Z);
    %surf(X,Y,abs(Zhat));
    %surf(X,Y,abs(fft2(img)));
    out = ifftshift(abs(ifft2(Zhat.* fft2(img))));
end

