clc
close all
clearvars

A=1; f_0=1;

f_s=200*f_0;

SNR=1/2;  %SNR=-3dB
f_c=f_s/10;     %cut-off frequency of the digital filter
[b,a]=butter(9,f_c/(f_s/2));  %designs a 9th-order low-pass digital Butterworth filter (IIR),where b is a vector containing coefficients of a moving average part and a is a vector containing coefficients of an auto-regressive part of the transfer function (see Equation (6.12) of Shin's book).

T_0=1/f_0;
T_vec=[2,100,1000]*T_0;
for T=T_vec
    [D_t,K]=samplingParameters_T_fs(T,f_s);
    t_vec=(0:K-1)*D_t;

    x_vec=A*sin(2*pi*f_0*t_vec);

    rng(0);
    n_vec=randn(1,K);   %create broad-band white noise
    n_vec=filtfilt(b,a,n_vec);  %performs zero-phase digital filtering. (Appendix H of Shin's book) The resulting sequence n is the band-limited (zero to f_c) white noise.
    n_vec=std(x_vec)/std(n_vec)/sqrt(SNR)*n_vec;
    x_hat_vec=x_vec+n_vec;

    r_xx_hat=xcorr(x_hat_vec,x_hat_vec,'unbiased');
    tau=(-(K-1):K-1)*D_t;

    figure
    subplot(2,1,1)
    plot(t_vec,x_hat_vec,t_vec,x_vec)
    xlim([0,2*T_0])
    xlabel('$t$ (sec.)','interpreter','latex');
    legend({'$\hat{x}(t)=x(t)+n(t)$','$x(t)$'},'interpreter','latex')
    title(['$\mathrm{SNR}=',num2str(SNR),',T=',num2str(T),'$ sec'],'interpreter','latex')

    %Additional optimization of the axes for correct comparison with the correlation curves
    pos=get(gca,'Position');
    pos(1)=pos(1)+pos(3)/2;
    pos(3)=pos(3)/2;
    set(gca,'Position',pos,'XGrid','on')

    subplot(2,1,2)
    plot(tau,r_xx_hat)
    xlim(2*T_0*[-1,1])
    xlabel('$\tau$ (sec.)','interpreter','latex');
    ylabel('$r_{\hat{x}\hat{x}}(\tau)$','interpreter','latex');
    set(gca,'XGrid','on')
end

export_figure(1:length(T_vec),'==',{'s1','s2','s3'})
