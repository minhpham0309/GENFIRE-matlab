a = normrnd(0,1,[500,500,500]);
tic;
b = fftshift(fftn((ifftshift(a))));
toc
clear b
%%
a = normrnd(0,1,[500,500,500]);
tic;
b = my_fft(a);
toc
clear b