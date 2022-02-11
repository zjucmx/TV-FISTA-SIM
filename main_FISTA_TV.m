close all;clear;clc;

%% reconstruction

load noiseimage.mat   % noiseimage is the raw data 
load OTF0.mat;  % OTF of 0 
load OTF1.mat;  % OTF of +-1 
load OTF2.mat;  % OTF of +-2
load I1.mat;    % lateral pattern of +-1
load I2.mat;    % lateral pattern of +-2
[xsize,ysize,zsize,r_num,p_num] = size(noiseimage);
% r_num is the number of different angles
% p_num is the number of different phases

wide = zeros(xsize,ysize,zsize);    % wide field
for i = 1:r_num
    for j =1:p_num
        wide = wide + noiseimage(:,:,:,i,j);
    end
end
wide = wide./max(wide(:));
xk = wide;
vk = xk;

para.maxiter = 100;         % iter number
para.tv_lambda =  0*1*0.003;    % TV weight 
para.tv_lambday =  1;       % TV weight in y (compared to x)
para.tv_lambdaz =  1;       % TV weight in z (compared to x)

lambda = para.tv_lambda;
prox = @(x)(1/2*(max(x,0)+tv3d(x, para.tv_lambda, para.tv_lambday, para.tv_lambdaz)));


% put everything in GPU
noiseimage = gpuArray(single(noiseimage));
I0 = 1;
I1 = gpuArray(single(I1));
I2 = gpuArray(single(I2));
OTF0 = gpuArray(single(OTF0));
OTF1 = gpuArray(single(OTF1));
OTF2 = gpuArray(single(OTF2));

xk = gpuArray(single(xk));
vk = gpuArray(single(vk));

%% 1. calculater L through power iteration
L_it  = zeros(r_num,p_num);
power_it = 6;

x = gpuArray(single(ones(xsize,ysize,zsize)));
for i = 1:r_num
    for j = 1:p_num
        
        for it = 1:power_it
            temp0 = x.* I0; temp0 = fftshift(fftn(temp0));temp0 = temp0.*OTF0;
            temp1 = x.* I1(:,:,:,i,j); temp1 = fftshift(fftn(temp1));temp1 = temp1.*OTF1;
            temp2 = x.* I2(:,:,:,i,j); temp2 = fftshift(fftn(temp2));temp2 = temp2.*OTF2;
            x = temp0+temp1+temp2;
        end
        
        temp0 = x.* I0; temp0 = fftshift(fftn(temp0));temp0 = temp0.*OTF0;
        temp1 = x.* I1(:,:,:,i,j); temp1 = fftshift(fftn(temp1));temp1 = temp1.*OTF1;
        temp2 = x.* I2(:,:,:,i,j); temp2 = fftshift(fftn(temp2));temp2 = temp2.*OTF2;
        x_up = temp0+temp1+temp2;
        temp_up = x_up.*x;
        temp_low = x.*x;
        L_it(i,j) = gather(sum(temp_up(:))/sum(temp_low(:)));
        L_it = real(L_it);
    end
end


%% 2. FISTA
for i = 1:r_num
    for j = 1:p_num
         
        
        L = L_it(i,j);
        tk = 1.0;
        for num=1:para.maxiter
                xold = xk;
                vold = vk;
                told = tk;

                imagef = fftshift(fftn(noiseimage(:,:,:,i,j)));
                temp0 = vk.* I0; temp0 = fftshift(fftn(temp0));temp0 = temp0.*OTF0;
                temp1 = vk.* I1(:,:,:,i,j); temp1 = fftshift(fftn(temp1));temp1 = temp1.*OTF1;
                temp2 = vk.* I2(:,:,:,i,j); temp2 = fftshift(fftn(temp2));temp2 = temp2.*OTF2;
                temp = temp0+temp1+temp2;
                temp = temp-imagef;

                k0 = (ifftn(ifftshift((temp.*OTF0)))).*I0;
                k1 = (ifftn(ifftshift((temp.*OTF1)))).*I1(:,:,:,i,j);
                k2 = (ifftn(ifftshift((temp.*OTF2)))).*I2(:,:,:,i,j);
                grads = 1/3*(k0+k1+k2);
                
                xk = prox(vold - 1/L*grads);
                xk(xk<0)=0; 
                tk = 1 + sqrt(1+4*told^2)/2;
                vk = xk + (told-1)/tk *(xk- xold);

        end
        
        
    end     
end
