clc
close all
clearvars

%Sampling parameters
f_s=1000;
N=500;
K=N;
[T,D_t,D_f]=samplingParameters_fs_N(f_s,N);
f_col=(0:N-1).'*D_f;

%Signal Parameters
A_signals=[1,2];
f_0_col=floor([.1;.3]*N+[0;.5])*D_f; %1st f_0 yields best case leakage while 2nd f_0 yields worst-case leakage
SNR=.5;

%Welch parapeters
alpha=0.5;
K_overlap=round(alpha*K)+1;
P=20;
K_tot=P*(K-K_overlap+1)+K_overlap-1;
t_tot_col=(0:K_tot-1).'*D_t;

%Generate signal
x_tot_col=(A_signals*sin(2*pi*f_0_col*t_tot_col.')).';
%Add noise to signal
rng(0);
x_tot_col=x_tot_col+std(x_tot_col)/sqrt(SNR)*randn(size(x_tot_col));

%Periodogram Estimator
R_XX=periodogram(x_tot_col(1:K),[],N,f_s,'twosided')*f_s;
semilogy(f_col,R_XX);
hold on

%Modified Periodogram Estimator
win_col=hann(K);
%win_col=ones(K,1);
R_XX=periodogram(x_tot_col(1:K),win_col,N,f_s,'twosided')*f_s;
semilogy(f_col,R_XX);

%Welch (averaged) cpsd estimation
R_XX_Welch_Matlab=cpsd(x_tot_col,x_tot_col,win_col,K_overlap,N,f_s,'twosided')*f_s;
semilogy(f_col,R_XX_Welch_Matlab)

xlim([0,f_s/2])
grid on
set(gca,'PlotBoxAspectRatio',[2,1,1])
xlabel('$f$ (Hz)','interpreter','latex');
legend_str={'$\hat{R}_{XX}^{\mathrm{periodo}}$','$\hat{R}_{XX}^{\mathrm{m-periodo}}$ using Hann window',['$\hat{R}_{XX}^{\mathrm{Welch}}$ for $P=',int2str(P),'$ using Hann window']};
legend(legend_str,'interpreter','latex')

export_figure(gcf,'',{'NonParametricSpectralEstimators'});