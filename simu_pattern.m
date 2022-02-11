close all;clear;clc;

r_num=3;
p_num=5;
lamda_ex=640;   %illumination wavelength  nm
lamda=670;  %fluorescence emission wavelength   nm
NA=1.49;
n2 = 1.515;  
pixelsize=38; 
k = 2*pi/lamda_ex*n2;
load OTF.mat  % OTF 
load sample.mat   % sample
sample = sample./max(sample(:));
rotation = [0 60 120]/180*pi;
phase_off = [0 0 0]/180*pi;


%% imaging model according to 4pi configuration
[xsize,ysize,zsize]=size(sample);
xc=floor(xsize/2+1); xx = -xc+1:1:-xc+xsize;
yc=floor(ysize/2+1); yy = -yc+1:1:-yc+ysize;
[xaxis,yaxis] = meshgrid(xx,yy);
xaxis = xaxis * pixelsize; yaxis = yaxis * pixelsize;
zc=floor(zsize/2+1); zz = -zc+1:1:-zc+zsize;
z_real = zz*pixelsize;

theta = pi/3; % the angle between adjacent beams
% Iz0 = 6+4*cos(2kcos(theta)z)+2cos(2kz)
% Iz1 = 4* cos( k(1-cos(theta))*z ) + 4*cos( k(1+cos(theta))  )
% Iz2 = 2+2*cos(2*k*cos(theta))
Iz0 = zeros(xsize,ysize,zsize);
Iz1 = zeros(xsize,ysize,zsize);
Iz2 = zeros(xsize,ysize,zsize);
for z = 1:1:zsize
    Iz0(:,:,z)=6+4*cos(2*k*cos(theta)*(z-zc)*pixelsize)+2*cos(2*k*(z-zc)*pixelsize);
    Iz1(:,:,z)=4*cos(k*(1-cos(theta))*(z-zc)*pixelsize)+4*cos(k*(1+cos(theta))*(z-zc)*pixelsize);
    Iz2(:,:,z)= 2+2*cos(2*k*cos(theta)*(z-zc)*pixelsize);
end

iPSF = abs(ifftn(ifftshift(OTF)));


H0 = iPSF.*Iz0;
H1 = iPSF.*Iz1;
H2 = iPSF.*Iz2;


OTF0 = abs(fftshift(fftn(H0)));
OTF1 = abs(fftshift(fftn(H1)));
OTF2 = abs(fftshift(fftn(H2)));
max0 = max(OTF0(:));
max1 = max(OTF1(:));
max2 = max(OTF2(:));
OTF0 = OTF0./max([max0,max1,max2]);
OTF1 = OTF1./max([max0,max1,max2]);
OTF2 = OTF2./max([max0,max1,max2]);


% calculater diffferent lateral pattern
I0 = 1;
I1 = zeros(xsize,ysize,zsize,r_num,p_num);
I2 = zeros(xsize,ysize,zsize,r_num,p_num);
for i =1:r_num
    for j = 1:p_num
        kx = k*sin(theta)*cos(rotation(i));
        ky = k*sin(theta)*sin(rotation(i));
        temp_I1 = cos(kx.*xaxis+ky.*yaxis+2*pi/p_num*(j-1)+phase_off(j));
        temp_I2 = cos(2*kx.*xaxis+2*ky.*yaxis+2*2*pi/p_num*(j-1)+2*phase_off(j));
        for slice = 1:zsize
            I1(:,:,slice,i,j) = temp_I1;
            I2(:,:,slice,i,j) = temp_I2;
        end
    end
end



noiseimage = zeros(xsize,ysize,zsize,r_num,p_num);
for i=1:r_num
    for j=1:p_num
        temp0 = sample*I0;temp0 = fftshift(fftn(temp0));temp0 = temp0.*OTF0;
        temp1 = sample.* I1(:,:,:,i,j);temp1 = fftshift(fftn(temp1));temp1 = temp1.*OTF1;
        temp2 = sample.* I2(:,:,:,i,j);temp2 = fftshift(fftn(temp2));temp2 = temp2.*OTF2;
        imagef = temp0 + temp1 + temp2;
        temp = abs(ifftn(ifftshift(imagef)));
        temp = temp + 0.05*randn(xsize,ysize,zsize)*max(temp(:));  % add different noise
        temp = temp - min(temp(:));
        noiseimage(:,:,:,i,j) = temp./max(temp(:));
    end
end


max0 = max(I0(:));
max1 = max(I1(:));
max2 = max(I2(:));
I1 = I1./max([max0,max1,max2]);
I2 = I2./max([max0,max1,max2]);


