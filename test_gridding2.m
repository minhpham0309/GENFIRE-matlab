img1 = imread('data\lena_128.jpg');
img1 = double(img1); img1 = img1/max(img1(:));

img2 = padarray(img1,[192,192],0,'both');

F = fftshift(fft2(img2));
FF = F(:);

[X,Y] = meshgrid(1:512,1:512);
XX = X(:);
YY = Y(:);
%%
XX_sample = [XX;XX] + 0.2*normrnd(0,1,[2*512^2,1]);
YY_sample = [YY;YY] + 0.2*normrnd(0,1,[2*512^2,1]);

F_sample = interp2(X,Y,F,XX_sample,YY_sample);
%%
F_grid = griddata(XX_sample,YY_sample,F_sample,XX,YY,'nearest');
%%
F_grid = reshape(F_grid,[512,512]);
figure; img(log(F_grid)); colorbar
F_grid(isnan(F_grid))=0;
img3 = ifft2(ifftshift(F_grid));
figure;img(img3);title('nearest');


