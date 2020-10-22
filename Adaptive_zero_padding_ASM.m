clc
clear all
close all

%% -----------------------------------------------
% This is a case of 2-D adaptive zero-padding ASM
%% -----------------------------------------------

%% input parameters
z = 8;                     % propagation distance/mm
n0 = 2000;                  % sampling number
pitch = 0.002;             % sampling interval/mm
lam = 500e-6;              % wavelength/mm
k = 2*pi/lam;
lb = n0*pitch;
xb = linspace(-lb/2,lb/2-lb/n0,n0)';

%% input field: rectangle aperture
t = zeros(n0,n0);
rr = 500;
t(n0/2-rr:n0/2+rr,n0/2-rr:n0/2+rr) = 1;
figure,imshow(t)

%% proposed adaptive zero padding ASM
N_p = round(lam*z/2/pitch^2/sqrt(1-(lam/2/pitch)^2));
if mod (N_p,2)==1
    N_p = N_p+1;
end

t_pad = padarray(t,[N_p/2,N_p/2]);
n = n0+N_p;
[fx,fy] = meshgrid(linspace(-1/2/pitch,1/2/pitch-1/n/pitch,n)',linspace(-1/2/pitch,1/2/pitch-1/n/pitch,n)');
H = exp(1i*k*z*sqrt(1-(lam*fx).^2-(lam*fy).^2));

tic
t_pad_FT = fftshift(fft2(fftshift(t_pad)));
t_pro = ifftshift(ifft2(ifftshift(t_pad_FT.*H)));
toc


L = n*pitch;
X = linspace(-L/2,L/2-L/n,n);
aa = t_pro(n/2-n0/2:n/2+n0/2-1,n/2-n0/2:n/2+n0/2-1);
figure,imshow(abs(aa),[])
title('proposed method: ampliude')

%% conventional ASM

N_c = n0;
tc_pad = padarray(t,[N_c/2,N_c/2]);

nc = n0+N_c;
[fxc,fyc] = meshgrid(linspace(-1/2/pitch,1/2/pitch-1/nc/pitch,nc),linspace(-1/2/pitch,1/2/pitch-1/nc/pitch,nc));
Hc = exp(1i*k*z*sqrt(1-(lam*fxc).^2-(lam*fyc).^2));

tic
tc_pad_FT = fftshift(fft2(fftshift(tc_pad)));
t_con = ifftshift(ifft2(ifftshift(tc_pad_FT.*Hc)));
toc

Lc = nc*pitch;
Xc = linspace(-Lc/2,Lc/2-Lc/nc,nc);
bb = t_con(nc/2-n0/2:nc/2+n0/2-1,nc/2-n0/2:nc/2+n0/2-1);
figure,imshow(abs(bb),[])
title('conventional method: ampliude')


corr2(abs(aa),abs(bb))
corr2(angle(aa),angle(bb))




