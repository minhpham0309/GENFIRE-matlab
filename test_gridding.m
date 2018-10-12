img1 = imread('data\lena_128.jpg');
img1 = double(img1); img1 = img1/max(img1(:));

img2 = padarray(img1,[192,192],0,'both');
supp = img2~=0;
F = fftshift(fft2(img2));
FF = F(:);

[X,Y] = meshgrid(-255.5:1:255.5,-255.5:1:255.5);
XX = X(:);
YY = Y(:);
%%
%  cos -sin   x
%  sin  cos   0
line1 = (1:512)'-256.5;
n_pj = 80;
theta = pi/n_pj;
XX_sample=[];
YY_sample=[];
for i=1:n_pj
    XX_sample_i = line1*cos((i-1)*theta);
    YY_sample_i = line1*sin((i-1)*theta);
    XX_sample = [XX_sample;XX_sample_i];
    YY_sample = [YY_sample;YY_sample_i];
end


%XX_sample = [XX;XX] + 0.1*normrnd(0,1,[2*512^2,1]);
%YY_sample = [YY;YY] + 0.1*normrnd(0,1,[2*512^2,1]);
%%
F_sample = interp2(X,Y,F,XX_sample,YY_sample);
%% interpolation by common methods
F_grid = griddata(XX_sample,YY_sample,F_sample,XX,YY,'natural');
F_grid(isnan(F_grid))=0;
%% interpolation by mean curvature
[~,~,F_grid] = mincurvi(XX_sample,YY_sample,F_sample,XX,YY);
%% interpolation by kriging
[XX,YY,F_grid] = kriging(XX_sample,YY_sample,F_sample,X,Y);
%% delete bad interpolation
mask = zeros(512^2,1,'single');
for i=1:512^2
    x_i = XX(i); y_i = YY(i);
    distance_matrix = sqrt((XX_sample - x_i).^2 + (YY_sample - y_i).^2);
    if find(distance_matrix<0.5), mask(i)=1;end
end
mask = reshape(mask,512,512);
mask = mask~=0;
figure; img(mask)
%% inverse distance weight
for i=1:512^2
    x_i = XX(i); y_i = YY(i);
    distance_matrix = sqrt((XX_sample - x_i).^2 + (YY_sample - y_i).^2);
    index_i = find(distance_matrix<1);
    F_grid(i) = sum(1./distance_matrix(index_i).*F_sample(index_i)) / sum(1./distance_matrix(index_i));
end
%% algorithm
obj = rand(512,512,'single');
u = obj;
%u = img2;

[XX1,YY1] = meshgrid(1:512,1:512);
x_cen = floor(512/2);
y_cen = floor(512/2);
kernel = (XX1-x_cen).^2 + (YY1-y_cen).^2;
sigma = 512;
kernel = exp(-kernel/sigma^2);

for k=1:500
    z = fftshift(fftn(u));
    err = sum(abs(z(mask) - F_grid(mask)))/sum(abs(F_grid(mask)));
    if mod(k,10)==0, fprintf('%d.error = %f\n',k,err);end
    z(mask) = F_grid(mask) ;%+ 0.5*z(mask);
    u_k = ifftn(ifftshift(z));
    obj = 1.3*u_k - .3*u;
    obj = max(0,real(obj)).*supp;
    %F_obj = fftshift(fftn(obj)) .* kernel; obj = real(ifftn(ifftshift(F_obj)));
    
    u = obj + .3*(u - u_k);
end
figure; img(obj,'colormap','gray'); title('natural');
%%
F_grid = reshape(F_grid,[512,512]);
F_grid(isnan(F_grid))=0;
figure;img3 = ifft2(F_grid);
img(img3);title('nearest');
%%
figure; img((abs((reshape(F_grid,512,512))))); colorbar
figure; img(((F)));colorbar