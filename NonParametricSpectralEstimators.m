clc
close all
clearvars

rng(0);

SNR=0.5;  %SNR_db=-3
f_c_by_f_Nyquist=[];

f_s=1000;
N=500;
K=N;
[T,D_t,D_f]=samplingParameters_fs_N(f_s,N);
f_col=(0:N-1).'*D_f;

%Signal Parameters
A_row=[1,2];
f_0_col=(floor([.1;.3]*N)+[0;.5])*D_f; %1st f_0 yields best case leakage while 2nd f_0 yields worst-case leakage

%Welch parapeters
alpha=0.5;
K_overlap=round(alpha*K)+1;
P=200;
K_tot=P*(K-K_overlap+1)+K_overlap-1;
t_tot_col=(0:K_tot-1).'*D_t;

%Generate signal
x_col=(A_row*sin(2*pi*f_0_col*t_tot_col.')).';
x_hat_col=addNoise(x_col,SNR,f_c_by_f_Nyquist);

%Periodogram Estimator
R_XX_periodogram=periodogram(x_hat_col(1:K),[],N,f_s,'twosided')*f_s;

%Modified Periodogram Estimator
win_col=hann(K);
%win_col=ones(K,1);
R_XX_m_periodogram=periodogram(x_hat_col(1:K),win_col,N,f_s,'twosided')*f_s;

%Welch (averaged) cpsd estimation
R_XX_Welch=cpsd(x_hat_col,x_hat_col,win_col,K_overlap,N,f_s,'twosided')*f_s;

semilogy(f_col/f_s,[R_XX_periodogram,R_XX_m_periodogram,R_XX_Welch])
xlim([0,.5])
grid on
set(gca,'PlotBoxAspectRatio',[2,1,1])
xlabel('$f/f_{\mathrm{s}}$','interpreter','latex');
legend_str=["$\hat{R}_{XX}^{\mathrm{raw}}$","$\hat{R}_{\bar{X}\bar{X}}^{\mathrm{raw}}$ using Hann window","$\hat{R}_{XX}^{\mathrm{Welch}}$ using Hann window \& $P="+P+'$'];
legend(legend_str,'interpreter','latex','Location','southeast')

% export_figure(gcf,'',"NonParametricSpectralEstimators");