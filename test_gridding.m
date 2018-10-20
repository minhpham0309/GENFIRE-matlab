%addpath('gridding')
addpath('nufftall-1.33')
img1 = imread('data\lena_128.jpg');
img1 = double(img1); img1 = img1/max(img1(:));

img2 = padarray(img1,[192,192],0,'both');
supp = padarray(ones(size(img1)),[192,192],0,'both');
supp = logical(supp);
F = fftshift(fft2(ifftshift(img2)));
FF = F(:);

[X,Y] = meshgrid(-256:1:255,-256:1:255);
XX = X(:);
YY = Y(:);
%%
%  cos -sin   x
%  sin  cos   0
line1 = (1:512)'-256.5;
n_pj = 40;
theta = pi/n_pj;
XX_sample=[];
YY_sample=[];
for i=1:n_pj
    XX_sample_i = line1*cos((i-1)*theta)-0.5;
    YY_sample_i = line1*sin((i-1)*theta)-0.5;
    XX_sample = [XX_sample;XX_sample_i];
    YY_sample = [YY_sample;YY_sample_i];
end


%XX_sample = [XX;XX] + 0.1*normrnd(0,1,[2*512^2,1]);
%YY_sample = [YY;YY] + 0.1*normrnd(0,1,[2*512^2,1]);
%%
F_sample = interp2(X,Y,F,XX_sample,YY_sample,'linear');

%% non uniform fft
%  dirft2d2(nj,xj,yj,iflag=1,ms,mt,fk);
XX_sample = mod(XX_sample,512)-256;
YY_sample = mod(YY_sample,512)-256;
%%
xj=[XX;XX_sample]; xj=(xj/max(abs(xj)))*pi;
yj=[YY;YY_sample]; yj=(yj/max(abs(yj)))*pi;
[min(xj),max(xj)]
FF = (F);
fk = [FF(:);F_sample(:)];
nj = length(xj);

for k=1:1
u = nufft2d1(nj,yj,xj,fk,1,1e-6,512,512);
u2 = ifftn(ifftshift(FF));
u = u;
[max(u(:)), max(u2(:))]
u = max(0,real(u)).*((supp))/pi^2;
%max(u(:))
%u = u/512^2;
figure; img(u)
%u = fftshift(u);
%u(u<0.1&fftshift(supp)) = max(u(:));
fk = nufft2d2(nj,yj,xj,-1,1e-6,512,512,u);
f1 = fk(1:512^2);
f2 = fk(512^2+1:end);
z2 = fftshift(fft2(ifftshift(u)));
[max(z2(:)),max(f1),max(FF(:))]
[max(f2),max(F_sample)]
%fk = fk/512^2*nj;
fk(512^2+1:end)=F_sample;
end
%u = fftshift(u);
%u = u(193:320,193:320);

figure; img(u,'colormap','gray'); %caxis([0.6 1.9253])
figure; img(log(reshape(fk(1:512^2),512,512)));colorbar
figure; img(log(F));colorbar
[min(u(:)),max(u(:))]

%% interpolation by common methods
F_grid = griddata(XX_sample,YY_sample,F_sample,XX,YY,'linear');
F_grid(isnan(F_grid))=0;
F_grid = reshape(F_grid,512,512);
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
%%
sum(abs(F_grid(mask) - F(mask)))/sum(abs(F(mask)))
sum(abs(abs(F_grid(mask)) - abs(F(mask))))/sum(abs(F(mask)))
%% inverse distance weight
for i=1:512^2
    x_i = XX(i); y_i = YY(i);
    distance_matrix = sqrt((XX_sample - x_i).^2 + (YY_sample - y_i).^2);
    index_i = find(distance_matrix<1);
    F_grid(i) = sum(1./distance_matrix(index_i).*F_sample(index_i)) / sum(1./distance_matrix(index_i));
end
%% algorithm
measuredK = ifftshift(F_grid);
measuredK_mask = ifftshift(mask);
supp_ifft = ifftshift(supp);
obj = rand(512,512,'single');
u = obj;
%u = img2;

[XX1,YY1] = meshgrid(1:512,1:512);
x_cen = floor(512/2);
y_cen = floor(512/2);
kernel = (XX1-x_cen).^2 + (YY1-y_cen).^2;
sigma = 800;
kernel = exp(-kernel/sigma^2);

for k=1:200
    z = (fftn(u));
    err = sum(abs(abs(z(measuredK_mask)) - abs(measuredK(measuredK_mask))))/sum(abs(measuredK(measuredK_mask)));
    if mod(k,10)==0, fprintf('%d.error = %f\n',k,err);end
    z(measuredK_mask) = measuredK(measuredK_mask) ;%+ 0.5*z(mask);
    u_k = ifftn((z));
    obj = 1*u_k - .0*u;
    obj = max(0,real(obj)).*supp_ifft;
    %F_obj = fftshift(fftn(obj)) .* kernel; obj = real(ifftn(ifftshift(F_obj)));
    
    u = obj + .3*(u - u_k);
end
obj = fftshift(obj);
figure; img(obj,'colormap','gray'); title('linear');
%%
F_grid = reshape(F_grid,[512,512]);
F_grid(isnan(F_grid))=0;
figure;img3 = ifft2(F_grid);
img(img3);title('nearest');
%%
figure; img((abs((reshape(F_grid,512,512)))));title('gridding'); colorbar
figure; img(((F)));colorbar; title('original');
