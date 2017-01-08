function Example_8_3()
clc
close all

%DFT parameters
f_s=200;
T=4;
[D_t,K]=samplingParameters_T_fs(T,f_s);
t_vec=(0:K-1)*D_t;

%s(t) parameters
T_burst=T/8;
K_burst=K*T_burst/T;
s_burst=sin(2*pi*10*t_vec(1:K_burst));

%Generate x(t)
a=1;b=.8; c=.75;
D_1=1;D_2=2.5;
k_1=D_1/D_t;k_2=D_2/D_t;

x_vec=zeros(1,K);
x_vec(1:K_burst)=a*s_burst;

y_vec=zeros(1,K);
y_vec(k_1+1:k_1+K_burst)=y_vec(k_1+1:k_1+K_burst)+b*s_burst;
y_vec(k_2+1:k_2+K_burst)=y_vec(k_2+1:k_2+K_burst)+c*s_burst;

SNR_vec=[inf,1/sqrt(2)];
N_SNR=length(SNR_vec);
rng(0);
for ii=1:N_SNR
    figure
    %Add noise to x(t) & y(t)
    x_vec=x_vec+std(x_vec)*randn(1,K)/SNR_vec(ii);
    y_vec=y_vec+std(y_vec)*randn(1,K)/SNR_vec(ii);

    r_xy=xcorr(y_vec,x_vec,'unbiased'); %linear cross-correlation
    r_yy=xcorr(y_vec,y_vec,'unbiased'); %linear auto-correlation
    tau=(-(K-1):K-1)*D_t;

    subplot(4,1,1)
    plot(t_vec,x_vec)
    xlabel('$t$ (sec.)', 'interpreter', 'latex');
    ylabel('$x(t)$', 'interpreter', 'latex')
    title(['$\mathrm{SNR} = ',num2str(SNR_vec(ii)),'$'], 'interpreter', 'latex')
    
    %Additional optimization of the axes for correct comparison with the correlation curves
    pos=get(gca,'Position');
    pos(1)=pos(1)+pos(3)/2;
    pos(3)=pos(3)/2;
    set(gca,'Position',pos,'XGrid','on')
    
    subplot(4,1,2)
    plot(t_vec,y_vec)
    xlabel('$t$ (sec.)', 'interpreter', 'latex');
    ylabel('$y(t)$', 'interpreter', 'latex')

    %Additional optimization of the axes for correct comparison with the correlation curves
    pos=get(gca,'Position');
    pos(1)=pos(1)+pos(3)/2;
    pos(3)=pos(3)/2;
    set(gca,'Position',pos,'XGrid','on')


    subplot(4,1,3)
    plot(tau,r_yy)
    axis tight
    xlabel('$\tau$ (sec.)', 'interpreter', 'latex');
    ylabel('$r_{yy}(\tau)$', 'interpreter', 'latex');
    set(gca,'XGrid','on')
    
    subplot(4,1,4)
    plot(tau,r_xy)
    axis tight
    xlabel('$\tau$ (sec.)', 'interpreter', 'latex');
    ylabel('$r_{xy}(\tau)$', 'interpreter', 'latex')
    set(gca,'XGrid','on')
end

export_figure(1:N_SNR,'',{'SNR1','SNR2'})