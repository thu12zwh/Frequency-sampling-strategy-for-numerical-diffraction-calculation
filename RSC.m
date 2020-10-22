clc
clear all
close all


%% input parameters
z = 80;                     % propagation distance/mm
n0 = 500;                  % sampling number
pitch = 0.002;             % sampling interval/mm
lam = 500e-6;              % wavelength/mm
k = 2*pi/lam;
lb = n0*pitch;
xb = linspace(-lb/2,lb/2-lb/n0,n0)';

%% input field
maxf = 200;
t = chirp(xb,0,max(xb),maxf);
figure,plot(xb,(t))
title('input field')
%% n0-point zero-padded input field
t_pad = padarray(t,[n0/2,0]);
n = length(t_pad);
t_pad_FT = fftshift(fft(fftshift(t_pad)));


%% option 2 RSC
L = n*pitch;
xn = linspace(-L/2,L/2-L/n,n)';
r = sqrt(xn.^2+z^2);
h1 = 1/2/pi*z./r.*(1./r-1i*k).*exp(1i*k*r)./r;
H1 = fftshift(fft(fftshift(h1)));
t_pro = ifftshift(ifft(ifftshift(t_pad_FT.*H1)));
% figure,plot(xn,abs(t_pro)/max(abs(t_pro)))
% figure,plot(xn,angle(t_pro))
% figure,plot(imag(H1))

aa = t_pro(n/2-n0/2:n/2+n0/2-1);
figure,plot(xb,abs(aa)/max(abs(aa)))
title('option (2) RSC amplitude')

%% option 1 RSC
r2 = sqrt(xb.^2+z^2);
h2 = 1/2/pi*z./r2.*(1./r2-1i*k).*exp(1i*k*r2)./r2;
h2 = padarray(h2,[n0/2,0]);
H2 = fftshift(fft(fftshift(h2)));
t_con = ifftshift(ifft(ifftshift(t_pad_FT.*H2)));
figure,plot(xn,abs(t_con)/max(abs(t_con)))
title('option (1) RSC amplitude')

%% RSI
X = xn;
uu = zeros(n,1);

tic
for j = 1:n
    fun = @(xn) 1/2/pi*z./sqrt((X(j)-xn).^2+z^2).*(1./sqrt((X(j)-xn).^2+z^2)...
               -1i*k).*exp(1i*k*sqrt((X(j)-xn).^2+z^2))./sqrt((X(j)-xn).^2+z^2).*chirp(xn,0,max(xb),maxf);
           uu(j,1) = integral(fun,min(xb),max(xb));
end
toc

uu = uu/max(abs(uu));
amplitude_rsi = abs(uu);
phase_rsi = angle(uu);
figure,plot(xn,amplitude_rsi)
title('RSI amplitude')
