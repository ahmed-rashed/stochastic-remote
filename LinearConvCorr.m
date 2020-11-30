clc
close all
clearvars

set(groot,'DefaultAxesColorOrder',[1,0,0])
set(groot,'DefaultAxesLineStyleOrder','-o|-x|-^|-s')

x_vec=[1,.75,.5,.25,.25];
K=length(x_vec);
y_vec=[1,1,1,1,1];
Delta_t=1;
t_vec=(0:K-1)*Delta_t;
t_lin_conv_vec=0:2*K-2;
tau_lin_corr_vec=-(K-1):K-1;

%plot x and y signals
for ii=1:2
    subplot(3,2,ii);
    plot(t_vec,x_vec)
    
    subplot(3,2,2+ii);
    plot(t_vec,y_vec)
end
subplot(3,2,1);xlabel('$\kappa$','interpreter','latex');
subplot(3,2,2);xlabel('$k$','interpreter','latex');
subplot(3,2,3);xlabel('$\kappa$','interpreter','latex');
subplot(3,2,4);xlabel('$k$','interpreter','latex');

subplot(3,2,1);ylabel('$x_{\kappa}$','interpreter','latex')
subplot(3,2,2);ylabel('$x_k$','interpreter','latex')
subplot(3,2,3);ylabel('$h_{\kappa}$','interpreter','latex')
subplot(3,2,4);ylabel('$y_{k}$','interpreter','latex')

%calculate and plot convolution
subplot(3,2,5)
lin_conv_vec=ifft(fft(x_vec,2*K-1).*fft(y_vec,2*K-1))*Delta_t;
lin_conv_vec1=conv(x_vec,y_vec)*Delta_t;
lin_conv_vec2=cconv(x_vec,y_vec,2*K-1)*Delta_t;
plot(t_lin_conv_vec,lin_conv_vec,t_lin_conv_vec,lin_conv_vec1,t_lin_conv_vec,lin_conv_vec2)
xlabel('$k$','interpreter','latex')
ylabel('$\left(x\overline{*}h\right)_{k}$','interpreter','latex')
title('Linear Convolution');
legend({'FFT','conv','cconv'})

%calculate and plot correlation
subplot(3,2,6)
lin_corr_vec=slow_xcorr(x_vec,y_vec)*Delta_t;
lin_corr_vec1=ifft(fft(y_vec,2*K-1).*conj(fft(x_vec,2*K-1)))*Delta_t;
lin_corr_vec2=fliplr(conv(x_vec,fliplr(y_vec)))*Delta_t;
lin_corr_vec3=xcorr(y_vec,x_vec)*Delta_t;
plot(tau_lin_corr_vec,lin_corr_vec,t_lin_conv_vec,lin_corr_vec1,tau_lin_corr_vec,lin_corr_vec2,tau_lin_corr_vec,lin_corr_vec3);
xlabel('$\kappa$','interpreter','latex')
ylabel('$r_{xy,\kappa}^{\mathrm{lin}}$','interpreter','latex')
title('Linear Correlation');
legend({'sum','FFT','cconv','xcorr'})

%Additional optimization of the axes for correct comparison with the correlation curves
for ii=1:4
    pos=get(subplot(3,2,ii),'Position');
    pos(1)=pos(1)+pos(3)/3;
    pos(3)=pos(3)*1/3;
    set(subplot(3,2,ii),'Position',pos)
    ylim([0,1.1*max([x_vec(:);y_vec(:)])])
end
pos=get(subplot(3,2,5),'Position');
pos(1)=pos(1)+pos(3)/3;
pos(3)=pos(3)*2/3;
set(subplot(3,2,5),'Position',pos)

export_figure(gcf,'==',{'LinearConvCorr'})

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultAxesLineStyleOrder','remove')