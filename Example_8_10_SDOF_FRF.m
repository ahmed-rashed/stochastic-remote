clearvars
close all
clc


fs=100;
N=250;
K=N;
[T,D_t,D_f]=samplingParameters_fs_N(fs,N);
t=(0:N-1)*D_t;
f_row=(0:N-1)*D_f;



A=100; zeta=0.03; f_n=10; w_n=2*pi*f_n; w_d=sqrt(1-zeta^2)*w_n;
h=(A/w_d)*exp(-zeta*w_n*t).*sin(w_d*t);

ScaleFactor=5; % 40 and 800
TT=ScaleFactor*T;
NN=ScaleFactor*N
rng(0);
x=2*randn(1,NN); 

y=conv(h,x);
y=y(1:NN);
% y=filter(h,1,x);

kappa_max=K;
[r_xx_lin, kappa]=xcorr(x,x,kappa_max,'unbiased');
[r_yy_lin, kappa]=xcorr(y,y,kappa_max,'unbiased');
[r_xy_lin, kappa]=xcorr(y,x,kappa_max,'unbiased');
tau=kappa*D_t;

N_pad=2*kappa_max;
T_pad=2*T;
D_f_pad=D_f/2;
f_pad=(0:N_pad-1)*D_f_pad;

Sxx=fft(r_xx_lin(1:N_pad));
Syy=fft(r_yy_lin(1:N_pad))/(fs^2);
Sxy=fft(r_xy_lin(1:N_pad))/fs;

H1=Sxy./Sxx;
H=fft(h,N_pad)/fs;

figure(1)
subplot(3,1,1)
plot(tau,r_xx_lin)
xlabel('$\tau$','interpreter', 'latex'); ylabel('$r_{xx}(\tau)$','interpreter', 'latex')

%figure(2) 
subplot(3,1,2)
plot(tau,r_yy_lin)
xlabel('$\tau$','interpreter', 'latex'); ylabel('$r_{yy}(\tau)$','interpreter', 'latex')

%figure(3) 
subplot(3,1,3)
plot(tau,r_xy_lin)
xlabel('$\tau$','interpreter', 'latex'); ylabel('$r_{xy}(\tau)$','interpreter', 'latex')

figure(4)
subplot(3,1,1:2)
plot(f_pad(1:N), 20*log10(abs(H(1:N))),f_pad(1:N), 20*log10(abs(H1(1:N))),'r');
xlabel('$f$ (Hz)','interpreter', 'latex');
ylabel('$|H(f)|$ (dB)','interpreter', 'latex')
legend({'$H(f)$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')

%figure(5)
subplot(3,1,3)
plot(f_pad(1:N), unwrap(angle(H(1:N))),f_pad(1:N), unwrap(angle(H1(1:N))),'r');
xlabel('$f$ (Hz)','interpreter', 'latex');
ylabel('$\angle H(f)$ (rad)','interpreter', 'latex')
legend({'$H(f)$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')

%% Comments 1 (use T=2000) %%

% figure(1)
% plot(t,4*h); hold on
% plot(tau(maxlag+1:end), r_xy_lin(maxlag+1:end), 'r:'); hold off
% xlabel('Time (s) and lag (\it\tau\rm)'); ylabel('Amplitude')
% 
% figure(2)
% Rhh=xcorr(h,h,kappa_max);
% plot(tau,4*Rhh); hold on
% plot(tau, r_yy_lin, 'r:'); hold off
% xlabel('Lag (\it\tau\rm)'); ylabel('Amplitude')

%% Comments 2 (use T=100) %%
Sxx_circ_Corr=cpsd(x,x,ones(size(x)),0,N_pad, fs, 'twosided');
Sxy_circ_Corr=cpsd(y/fs,x, ones(size(y)),0, N_pad, fs, 'twosided');
H1_circ_Corr=Sxy_circ_Corr./Sxx_circ_Corr;

%figure(1)
figure(5)
subplot(3,1,1:2)
plot(f_pad(1:N), 20*log10(abs(H1_circ_Corr(1:N))),f_pad(1:N), 20*log10(abs(H1(1:N))), 'r');
xlabel('$f$ (Hz)','interpreter', 'latex');
ylabel('$|H(f)|$ (dB)','interpreter', 'latex')
legend({'$H_{1}^{\mathrm{circ corr}}$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')

%figure(2)
subplot(3,1,3)
plot(f_pad(1:N), unwrap(angle(H1_circ_Corr(1:N))),f_pad(1:N), unwrap(angle(H1(1:N))), 'r'); hold on;set(gca,'ColorOrderIndex',1);
xlabel('$f$ (Hz)','interpreter', 'latex');
ylabel('$\angle H(f)$ (rad)','interpreter', 'latex')
legend({'$H_{1}^{\mathrm{circ corr}}$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')

Sxx_w=cpsd(x,x, hanning(N_pad),N_pad/2, N_pad, fs, 'twosided');
Sxy_w=cpsd(y/fs,x, hanning(N_pad),N_pad/2, N_pad, fs, 'twosided');
H1_w=Sxy_w./Sxx_w;

figure(7)
subplot(3,1,1:2)
plot(f_pad(1:N), 20*log10(abs(H1_w(1:N))),f_pad(1:N), 20*log10(abs(H1(1:N))), 'r');
xlabel('$f$ (Hz)','interpreter', 'latex');
ylabel('$|H(f)|$ (dB)','interpreter', 'latex')
legend({'$H_{1}^{\mathrm{Welch}}$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')

%figure(2)
subplot(3,1,3)
plot(f_pad(1:N), unwrap(angle(H1_w(1:N))),f_pad(1:N), unwrap(angle(H1(1:N))), 'r'); hold on;set(gca,'ColorOrderIndex',1);
xlabel('$f$ (Hz)','interpreter', 'latex');
ylabel('$\angle H(f)$ (rad)','interpreter', 'latex')
legend({'$H_{1}^{\mathrm{Welch}}$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')
