function SpectralDensity()
close all
clc

%%%
f_s=1000;
D_t=1/f_s;

N=200;
T=N*D_t;

N_overlapRatio=0.5;
N_overlap=ceil(N*N_overlapRatio);
Q=20;
N_total=Q*(N-N_overlap)+N_overlap;

f_0=200;
t_row_total=(0:N_total-1)*D_t;
x=sin(2*pi*f_0*t_row_total);
y=x*2;
x_noisy=x+rand(1,N_total);
y_noisy=y+rand(1,N_total);


%Raw cpsd estimation
R_XX_noisy_raw=cpsd(x_noisy,x_noisy,rectwin(N_total),0,N_total,f_s,'twosided')*f_s;
R_YY_noisy_raw=cpsd(y_noisy,y_noisy,rectwin(N_total),0,N_total,f_s,'twosided')*f_s;
R_XY_noisy_raw=cpsd(y_noisy,x_noisy,rectwin(N_total),0,N_total,f_s,'twosided')*f_s;
D_f_total=f_s/N_total;
f_col_total=(0:N_total-1).'*D_f_total;
figure
subplot(3,2,1)
semilogy(f_col_total(1:floor(N_total/2)),R_XX_noisy_raw(1:floor(N_total/2)))
grid
title(['$\hat{R}_{XX}$ for $N=',int2str(N_total),'$'], 'interpreter', 'latex')
subplot(3,2,3)
semilogy(f_col_total(1:floor(N_total/2)),R_YY_noisy_raw(1:floor(N_total/2)))
grid
title(['$\hat{R}_{YY}$ for $N=',int2str(N_total),'$'], 'interpreter', 'latex')
subplot(3,2,5)
semilogy(f_col_total(1:floor(N_total/2)),abs(R_XY_noisy_raw(1:floor(N_total/2))))
grid
title(['$\hat{R}_{XY}$ for $N=',int2str(N_total),'$'], 'interpreter', 'latex')
xlabel('$f$ (Hz)', 'interpreter', 'latex')

%Welch (averaged) cpsd estimation
R_XX_noisy_Welc=cpsd(x_noisy,x_noisy,hann(N),N_overlap,N,f_s,'twosided')*f_s;
R_YY_noisy_Welc=cpsd(y_noisy,y_noisy,hann(N),N_overlap,N,f_s,'twosided')*f_s;
R_XY_noisy_Welc=cpsd(y_noisy,x_noisy,hann(N),N_overlap,N,f_s,'twosided')*f_s;
D_f=f_s/N;
f_col=(0:N-1).'*D_f;
subplot(3,2,2)
semilogy(f_col(1:floor(N/2)),R_XX_noisy_Welc(1:floor(N/2)))
grid
title(['$R_{XX}^{\textrm{Welch}}$ for $N=',int2str(N),'$ and $Q=',int2str(Q),'$'], 'interpreter', 'latex')
subplot(3,2,4)
semilogy(f_col(1:floor(N/2)),R_YY_noisy_Welc(1:floor(N/2)))
grid
title(['$R_{YY}^{\textrm{Welch}}$ for $N=',int2str(N),'$ and $Q=',int2str(Q),'$'], 'interpreter', 'latex')
subplot(3,2,6)
semilogy(f_col(1:floor(N/2)),abs(R_XY_noisy_Welc(1:floor(N/2))))
grid
title(['$\left|R_{XY}^{\textrm{Welch}}\right|$ for $N=',int2str(N),'$ and $Q=',int2str(Q),'$'], 'interpreter', 'latex')
xlabel('$f$ (Hz)', 'interpreter', 'latex')

export_figure(gcf,'||',{'SpectralDensity'})

figure
H_XY_raw=tfestimate(x_noisy,y_noisy,rectwin(N_total),0,N_total,f_s,'twosided');
gamma_2_XY_raw=mscohere(x_noisy,y_noisy,rectwin(N_total),0,N_total,f_s,'twosided');
subplot(2,2,1)
semilogy(f_col_total(1:floor(N_total/2)),abs(H_XY_raw(1:floor(N_total/2))))
grid
title(['$\left|H_{YX}\right|$ for $N=',int2str(N_total),'$'], 'interpreter', 'latex')
subplot(2,2,3)
plot(f_col_total(1:floor(N_total/2)),gamma_2_XY_raw(1:floor(N_total/2)))
ylim([0,1.1]);grid
title(['$\gamma_{YX}^2$ for $N=',int2str(N_total),'$'], 'interpreter', 'latex')
xlabel('$f$ (Hz)', 'interpreter', 'latex')

H_XY_Welc=tfestimate(x_noisy,y_noisy,hann(N),N_overlap,N,f_s,'twosided');
gamma_2_XY_Welc=mscohere(x_noisy,y_noisy,hann(N),N_overlap,N,f_s,'twosided');
subplot(2,2,2)
semilogy(f_col(1:floor(N/2)),abs(H_XY_Welc(1:floor(N/2))))
grid
title(['$\left|H_{YX}\right|$ for $N=',int2str(N),'$ and $Q=',int2str(Q),'$'], 'interpreter', 'latex')
subplot(2,2,4)
plot(f_col(1:floor(N/2)),gamma_2_XY_Welc(1:floor(N/2)))
ylim([0,1.1]);grid
title(['$\gamma_{YX}^2$ for $N=',int2str(N),'$ and $Q=',int2str(Q),'$'], 'interpreter', 'latex')
xlabel('$f$ (Hz)', 'interpreter', 'latex')

export_figure(gcf,'',{'coherence'})