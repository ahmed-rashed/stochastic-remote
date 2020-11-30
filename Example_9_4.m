function Example_9_4()
clc
close all

set(groot,'DefaultAxesColorOrder',[0,0,1;1,0,1;1,0,0;0,0.5,0;0,0,0])
set(groot,'DefaultLineLineWidth',1);

A_r=[-0.31871,-0.15923];
f_r=[5,15];
zeta_r=[0.05,0.03];

f_s=50/15*max(f_r);
T=10;
[D_t,K,D_f]=samplingParameters_T_fs(T,f_s);
N=K;
t_col=(0:K-1).'*D_t;
f_col=(0:N-1).'*D_f;

K_Tot=5e3*K;   %T_Tot=5e3*T;  %Measured very long samples to do alot of averaging to minimize random errors on the estimation of the spectral density functions. This will be discussed in Chapter 10.


w_r=2*pi*f_r;
w_d_r=sqrt(1-zeta_r.^2).*w_r;
h_exact=(-2*A_r*(exp(-(zeta_r.*w_r).'*t_col.').*sin(w_d_r.'*t_col.'))).';
H_exact=(-2*(w_d_r.*A_r)*(1./(w_r.'.^2*ones(1,N)-ones(2,1)*(2*pi*f_col.').^2+2*1i*(zeta_r.*w_r).'*2*pi*f_col.'))).';
H_DFT=fft(h_exact)/f_s;

SNR_x=[2,inf,2];
SNR_y=[inf,4,4];
figs=nan(size(SNR_x));
rng(0);
x_vec_long=randn(1,K_Tot);
% y_vec_long=filter(h_exact,1,x_vec_long)/f_s; %In this case kappa(f)~=1 and hence equation (9.68) is not valid
y_vec_long=filter(h_exact,1,x_vec_long); %Not scaled to show the advantage of H_T. In this case kappa(f)=1 and hence equation (9.68) is valid
std_x=std(x_vec_long);
std_y=std(y_vec_long);
for nn=1:length(SNR_x)
    rng(10);
    x_noisy_vec_long=x_vec_long+std_x*randn(1,K_Tot)/SNR_x(nn);

    rng(20);
    y_noisy_vec_long=y_vec_long+std_y*randn(1,K_Tot)/SNR_y(nn);

    R_XX=cpsd(x_noisy_vec_long,x_noisy_vec_long,hanning(K),K/2,K,f_s,'twosided')*f_s; %Welch with 50% overlap
    R_YY=cpsd(y_noisy_vec_long,y_noisy_vec_long,hanning(K),K/2,K,f_s,'twosided')*f_s; %Welch with 50% overlap
    R_XY=cpsd(y_noisy_vec_long,x_noisy_vec_long,hanning(K),K/2,K,f_s,'twosided')*f_s; %Welch with 50% overlap
    R_YX=conj(R_XY);

    figs(nn)=figure;
    gamma_YX_2=abs(R_XY).^2./R_XX./R_YY;
    subplot(6,1,1:2);
    plot(f_col,gamma_YX_2)
    axis([0,f_s/2,0,1.1])
    grid
    set(gca,'XTickLabel',[]);
    ylabel('$\gamma_{YX}^{2}(f)$','interpreter','latex')
    title(['SNR_x=',num2str(SNR_x(nn)),' & SNR_y=',num2str(SNR_y(nn))])

    H_1=R_XY./R_XX; 
    H_2=R_YY./R_YX;
    H_T=(R_YY-R_XX + sqrt((R_XX-R_YY).^2 + 4*abs(R_XY).^2))./(2*R_YX);  %equation (9.68)

    ax_mag_h=subplot(6,1,3:5);
    ax_phase_h=subplot(6,1,6);
    plot_FRF_mag_phase(f_col,[H_exact,H_DFT,H_1,H_2,H_T],false,ax_mag_h,ax_phase_h);
    legend(ax_mag_h,{'$H_{\mathrm{exact}}$','$\mathrm{DFT}\left[h_{\mathrm{exact}}\right]$','$H_{1}$','$H_{2}$','$H_{\mathrm{T}}$'},'interpreter','latex')
    xlim(ax_mag_h,[0,f_s/2])
    xlim(ax_phase_h,[0,f_s/2])

    %Comments 2:
    [b,a]=invfreqz(H_T(1:floor(N/2)+1),2*pi*f_col(1:floor(N/2)+1)/f_s,4,4,[],30); 
    Hz=freqz(b,a,floor(N/2)+1,f_s);
    figure
    plot_FRF_mag_phase(f_col(1:floor(N/2)+1),Hz,false,[],[],'','Hz');
end

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultLineLineWidth','remove');

export_figure(figs,'||',{'H_estimators_1','H_estimators_2','H_estimators_3'})
