function Example_8_5()
clc
close all

set(groot,'DefaultAxesColorOrder',[0,0,1;1,0,0;0,0.5,0;1,0,1;0,0,0])
set(groot,'DefaultLineLineWidth',1);

A=1; f_0=1; T_0=1/f_0;
f_s=200*f_0;
% T=100*T_0;
T_vec=[2,100,1000]*T_0;
for T=T_vec
    [D_t,K]=samplingParameters_T_fs(T,f_s);
    t_vec=(0:K-1)*D_t;

    s_vec=A*sin(2*pi*f_0*t_vec);

    rng(0);
    SNR=1/sqrt(2);
    n=randn(1,K);
    f_c=f_s/10;
    [b,a] = butter(9,f_c/(f_s/2));
    n=filtfilt(b,a,n);
    n=std(s_vec)*n/std(n)/SNR; %SNR = -3dB
    y=s_vec+n;

    r_yy=xcorr(y,y,'unbiased');
    tau=(-(K-1):K-1)*D_t;

    figure
    subplot(2,1,1)
    plot(t_vec,y,t_vec,s_vec)
    xlim([0,2*T_0])
    xlabel('$t$ (sec.)', 'interpreter', 'latex');
    legend({'$y(t)=s(t)+n(t)$','$s(t)$'}, 'interpreter', 'latex')
    title(['$\textrm{SNR} = ',num2str(SNR),', T=',num2str(T),'$ sec'], 'interpreter', 'latex')

    %Additional optimization of the axes for correct comparison with the correlation curves
    pos=get(gca,'Position');
    pos(1)=pos(1)+pos(3)/2;
    pos(3)=pos(3)/2;
    set(gca,'Position',pos,'XGrid','on')

    subplot(2,1,2)
    plot(tau,r_yy)
    xlim(2*T_0*[-1,1])
    xlabel('$\tau$ (sec.)', 'interpreter', 'latex');
    ylabel('$r_{yy}(\tau)$', 'interpreter', 'latex');
    set(gca,'XGrid','on')
end

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultLineLineWidth','remove');

export_figure(1:length(T_vec),'',{'s1','s2','s3'})