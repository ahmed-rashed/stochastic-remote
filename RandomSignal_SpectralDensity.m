clc
close all
clearvars

%DFT parameters
a=1; b=0.8; c=0.75;
D_T_1=1;
D_T_2=3*D_T_1;
D_T=abs(D_T_2-D_T_1);
SNR_y=10;

%We want to calculate r_{xy}(\tau) to up to \tau=2*max(D_T_1,D_T_2)
%Thus, 2*max(D_T_1,D_T_2) should be <= 0.1*T_Tot
%or D_T_2<= 0.05*T_Tot or T_Tot>=20*D_T_2
T_Tot=720*D_T_2;
f_s=20*2/D_T;   %To detect the delayed signals, f_s>2/D_T
% T_Tot=720*D_T_2;
% f_s=20*2/D_T;   %To detect the delayed signals, f_s>2/D_T
[D_t,K_Tot,D_f]=samplingParameters_T_fs(T_Tot,f_s);
f_c=0.4*f_s;
if f_c<(1/D_T),warning('s(t) is not wide band enough for detecting D_T=D_T_1 & D_T_2'),end

%s(t) parameters
T_start=max(D_T_1,D_T_2);
T_s=T_Tot+T_start;
K_s=T_s/D_t;    %K_Tot+k_start
rng(0);
s_vec=randn(1,K_s);   %white noise
[b_filt,a_filt]=butter(9,f_c/(f_s/2));  %designs a 9th-order low-pass digital Butterworth filter (IIR), where b is a vector containing coefficients of a moving average part and a is a vector containing coefficients of an auto-regressive part of the transfer function (see Equation (6.12) of Shin's book).
s_vec=filtfilt(b_filt,a_filt,s_vec);  %Filter s(t) with (full) bandwidth approximately 20 Hz (- fc to fc).
s_vec=s_vec-mean(s_vec);    % Makes mean(s)=0
s_vec=s_vec/std(s_vec);     % Makes std(s)=1

k_start=T_start/D_t;
x_vec=a*s_vec((1:K_Tot)+k_start);
k_1=D_T_1/D_t;k_2=D_T_2/D_t;
y_vec=b*s_vec((1:K_Tot)+k_start-k_1)+c*s_vec((1:K_Tot)+k_start-k_2);

rng(10);
y_vec=y_vec+randn(1,K_Tot)/SNR_y;

r_xy=xcorr(y_vec,x_vec,'unbiased');
tau=(-(K_Tot-1):K_Tot-1)*D_t;

figure
plot(tau,r_xy)
set(gca,'XGrid','on')
xlabel('$\tau$ (sec.)','interpreter','latex')
ylabel('$r_{xy}(\tau)$','interpreter','latex')
xlim([0,2*max(D_T_1,D_T_2)])

%Spectral Estimation
K=K_Tot/10;
N=K;
[~,~,D_f]=samplingParameters_fs_N(f_s,K);
%Welch with hanning window and 50% overlap
R_XX=cpsd(x_vec,x_vec,hanning(K),K/2,K,f_s,'twosided')*f_s;
R_YY=cpsd(y_vec,y_vec,hanning(K),K/2,K,f_s,'twosided')*f_s;
R_XY=cpsd(y_vec,x_vec,hanning(K),K/2,K,f_s,'twosided')*f_s;
f_vec=(0:N-1).'*D_f;

figure
ax1=subplot(4,1,1);
semilogy(f_vec(1:floor(N/2)),R_XX(1:floor(N/2)));
ylabel('$\hat{R}_{XX}(f)$','interpreter','latex')
set(ax1,'XGrid','on','XTickLabel',[],'YLimSpec','tight');
ax2=subplot(4,1,2);
semilogy(f_vec(1:floor(N/2)),R_YY(1:floor(N/2)));
ylabel('$\hat{R}_{YY}(f)$','interpreter','latex')
set(ax2,'XGrid','on','XTickLabel',[],'YLimSpec','tight');

ax_mag_h=subplot(4,1,3);
ax_phase_h=subplot(4,1,4);
plot_FRF_mag_phase(f_vec(1:floor(N/2)),R_XY(1:floor(N/2)),false,ax_mag_h,ax_phase_h,'','\hat{R}_{XY}');
set(ax_mag_h,'YLimSpec','tight')
set(ax_phase_h,'YLimSpec','tight')

export_figure(gcf,'==',{'Echo_Spectrum'})