function Example_8_10()

% Although the rectangular window is included in Table 10.2, it is rarely used since its spectral
% window side lobes cause large leakage effects.However, in Chapter 8, the rectangular window
% function is applied to estimate the spectral density functions in MATLAB Examples 8.8-8.10
% (note that we have used a very long data length T).

clc
close all

A=100; zeta=0.03; f_n=10; w_n=2*pi*f_n; w_d=sqrt(1-zeta^2)*w_n;

f_s=10*f_n;
N=250;
K=N;
D_f=f_s/N;
D_t=1/f_s;
T=1/D_f;
t_row=(0:K-1)*D_t;
f_row=(0:N-1)*D_f;

h_exact=A*exp(-zeta*w_n*t_row).*sin(w_d*t_row)/w_d;
H_exact=A./(-(2*pi*f_row).^2+2*1i*zeta*w_n*(2*pi*f_row)+w_n^2);
%H_fft=fft(h)/N*2;
H_fft=fft(h_exact)/f_s;

figure
[ax_mag_h,ax_phase_h]=plot_FRF_mag_phase(f_row(1:floor(N/2)), [H_exact(1:floor(N/2));H_fft(1:floor(N/2))],0);
legend(ax_mag_h,{'$H(f)^{\mathrm{exact}}$','$H^{\mathrm{FFT}}$'},'interpreter','latex')

RecordLengthMultiplier=800; % 40 and 800
TT=RecordLengthMultiplier*T;
KK=RecordLengthMultiplier*K;
rng(0);
x=2*randn(1,KK);    %std(x)=2

y=conv(h_exact,x);    %Linear convolution. This yields y of length KK+K+1
y=y(1:KK)*D_t;  %Use same number of elements as x
% y=filter(h,1,x);

[r_xx_lin, kappa]=xcorr(x,x,K,'unbiased');    %r_xx_lin is 2KK+1 length
[r_yy_lin, kappa]=xcorr(y,y,K,'unbiased');    %r_yy_lin is 2KK+1 length
[r_xy_lin, kappa]=xcorr(y,x,K,'unbiased');    %r_xy_lin is 2KK+1 length
tau=kappa*D_t;

figure
subplot(3,1,1)
plot(tau,r_xx_lin)
title('$r_{xx}^{\mathrm{lin corr}} (\tau)$','interpreter', 'latex')

subplot(3,1,2)
plot(tau,r_yy_lin)
title('$r_{yy}^{\mathrm{lin corr}} (\tau)$','interpreter', 'latex')

subplot(3,1,3)
plot(tau,r_xy_lin)
title('$r_{xy}^{\mathrm{lin corr}} (\tau)$','interpreter', 'latex')
xlabel('$\tau$','interpreter', 'latex');

K_pad=2*K-1;
N_pad=2*N-1;
T_pad=K_pad/K*T;
df_pad=N/N_pad*D_f;
f_pad=(0:N_pad-1)*df_pad;

R_XX_lin=fft(r_xx_lin(1:end-1))/(2*K-1)*2;    %R_XX_lin is 2K length
R_YY_lin=fft(r_yy_lin(1:end-1))/(2*K-1)*2;    %R_YY_lin is 2K length
R_XY_lin=fft(r_xy_lin(1:end-1))/(2*K-1)*2;    %R_XY_lin is 2K length

H1=R_XY_lin./R_XX_lin;
H_interpolated=fft(h_exact,N_pad)/N*2;
H_exact1=A./(-(2*pi*f_pad).^2+2*1i*zeta*w_n*(2*pi*f_pad)+w_n^2);

figure
[ax_mag_h,ax_phase_h]=plot_FRF_mag_phase(f_pad(1:N), [H_exact1(1:N);H_interpolated(1:N);H1(1:N)],0);
legend(ax_mag_h,{'$H(f)^{\mathrm{exact}}$','$H(f)^{\mathrm{interpolated}}$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')

%%Comments 1 (use T=2000) %%
figure
plot(t_row,h_exact,tau, r_xy_lin);
legend({'$h(t)$','$r_{xy}(\tau)$'},'interpreter','latex')
xlabel('$t$ \& $\tau$ (s)','interpreter','latex');
ylabel('Amplitude')

% figure(2)
% Rhh=xcorr(h,h,N);
% plot(tau,4*Rhh); hold on
% plot(tau, Ryy, 'r:'); hold off
% xlabel('Lag (\it\tau\rm)'); ylabel('Amplitude')



% Periodogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [R_XX_circ_Corr,ff]=cpsd(x,x,rectwin(KK),0,KK, f_s, 'twosided');
% [R_XY_circ_Corr,ff]=cpsd(y,x,rectwin(KK),0,KK, f_s, 'twosided');
% H1_circ_Corr=R_XY_circ_Corr./R_XX_circ_Corr;
[R_XX_circ_Corr,ff]=cpsd(x(1:K),x(1:K),rectwin(K),0,K, f_s, 'twosided');
[R_XY_circ_Corr,ff]=cpsd(y(1:K),x(1:K),rectwin(K),0,K, f_s, 'twosided');
H1_circ_Corr=R_XY_circ_Corr./R_XX_circ_Corr;

% figure
% [ax_mag_h,ax_phase_h]=plot_FRF_mag_phase(ff(1:floor(KK/2)), H1_circ_Corr(1:floor(KK/2)),0);
% legend(ax_mag_h,{'$H_{1}^{\mathrm{circ corr}}$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')
% grid(ax_mag_h,'off')
% grid(ax_phase_h,'off')

% Welch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_XX_Welch=cpsd(x,x, hanning(N_pad),floor(N_pad/2), N_pad, f_s, 'twosided').';
R_XY_Welch=cpsd(y,x, hanning(N_pad),floor(N_pad/2), N_pad, f_s, 'twosided').';
H1_Welch=R_XY_Welch./R_XX_Welch;

figure
[ax_mag_h,ax_phase_h]=plot_FRF_mag_phase(f_pad(1:N), [H1_circ_Corr.';H1_Welch(1:N);H1(1:N)],0);
legend(ax_mag_h,{'$H_{1}^{\mathrm{periodogram}}$','$H_{1}^{\mathrm{Welch}}$','$H_{1}^{\mathrm{lin corr}}$'},'interpreter','latex')
grid(ax_mag_h,'off')
grid(ax_phase_h,'off')
