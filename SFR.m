clc
clear all
close all

%% input parameters


z = 10;                    % propagation distance/mm
n0 = 500;                  % sampling number
pitch = 0.002;             % sampling interval/mm
lam = 500e-6;              % wavelength/mm
k = 2*pi/lam;
lb = n0*pitch;
xb = linspace(-lb/2,lb/2-lb/n0,n0)';

f_max = 1/2/pitch-lb/2/lam/z;  % maximum spatial frequency without aliasing (Eq. (6))

%% input field: chirp grating
maxf = 200;
t = chirp(xb,0,max(xb),maxf);
figure,plot(xb,(t))
title('input field')

%% sampling number of output field N_hat
nn = floor(lam*z/pitch^2-n0);
if nn<n0
    nn = n0;
end


%% conventional method
l0 = n0*pitch;
x0 = linspace(-l0/2,l0/2-l0/n0,n0)';
chirp0 = exp(1i*k/2/z*x0.^2);
t_c0 = t.*chirp0;
t_pro0 = fftshift(fft(fftshift(padarray(t_c0,0))));
t_pro0 = t_pro0./max(abs(t_pro0));
L0 = lam*z/pitch;
X0 = linspace(-L0/2,L0/2-L0/n0,n0)';
Chirp0 = exp(1i*k/2/z*X0.^2);
t_con = Chirp0.*t_pro0;

% full range result
figure,plot(X0,abs(t_con))
title('conventional method amplitude')
figure,plot(X0,angle(t_con))
title('conventional method phase')

% central range result(Eq. 8)
LL = lam*z/pitch-n0*pitch;
aa = t_con(abs(X0)<=LL/2);
n1 = length(aa);
xx = linspace(-LL/2,LL/2,n1);
figure,plot(xx,abs(aa))
title('conventional method amplitude:central part')
figure,plot(xx,angle(aa))
title('conventional method phase:central part')

%% proposed method

t_pron = fftshift(fft(fftshift(padarray(t_c0,[(nn-n0)/2,0]))));
t_pron = t_pron./max(abs(t_pron));
Xn = linspace(-L0/2,L0/2-L0/nn,nn)';
Chirpn = exp(1i*k/2/z*Xn.^2);
t_pro = Chirpn.*t_pron;

bb = t_pro(abs(Xn)<=LL/2);
n2 = length(bb);
xxx = linspace(-LL/2,LL/2,n2);
figure,plot(xxx,abs(bb))
title('proposed method amplitude')
figure,plot(xxx,angle(bb))
title('proposed method phase')


%% analytical Fresnel integral 
L = lam*z/pitch;
X = linspace(-L/2,L/2-L/nn,nn)';
uu = zeros(nn,1);

tic
for j = 1:nn
    fun = @(xn) exp(1i*k/2/z*(X(j)-xn).^2).*chirp(xn,0,max(xb),maxf);
           uu(j,1) = integral(fun,min(xb),max(xb));
end
toc

uu = uu/max(abs(uu));
cc = uu(abs(X)<=LL/2);
n3 = length(cc);
xxxx = linspace(-LL/2,LL/2,n3);
figure,plot(xxxx,abs(cc))
title('Analytical integral amplitude')
figure,plot(xxxx,angle(cc))
title('Analytical integral phase')

%% correlation
corr(abs(bb),abs(cc))
corr((angle(bb)),(angle(cc)))
