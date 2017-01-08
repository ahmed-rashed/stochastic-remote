clc
close all
clearvars

x_vec=[1,.75,.5,.25,.25];
K=length(x_vec);
y_vec=[1,1,1,1,1];
t_vec=0:K-1;
t_lin_conv_vec=0:2*K-2;
tau_lin_corr_vec=-(K-1):K-1;

%plot x and y signals
for ii=1:2
    subplot(3,2,ii);
    plot(t_vec,x_vec,'-o')
    
    subplot(3,2,2+ii);
    plot(t_vec,y_vec,'-o')
end
subplot(3,2,1);xlabel('$\kappa$', 'interpreter', 'latex');
subplot(3,2,2);xlabel('$k$', 'interpreter', 'latex');
subplot(3,2,3);xlabel('$\kappa$', 'interpreter', 'latex');
subplot(3,2,4);xlabel('$k$', 'interpreter', 'latex');

subplot(3,2,1);ylabel('$x_{\kappa}$', 'interpreter', 'latex')
subplot(3,2,2);ylabel('$x_k$', 'interpreter', 'latex')
subplot(3,2,3);ylabel('$h_{\kappa}$', 'interpreter', 'latex')
subplot(3,2,4);ylabel('$y_{k}$', 'interpreter', 'latex')

%calculate and plot convolution
subplot(3,2,5)
lin_conv_vec=ifft(fft(x_vec,2*K-1).*fft(y_vec,2*K-1));plot(t_lin_conv_vec,lin_conv_vec,'-o');hold on
lin_conv_vec1=conv(x_vec,y_vec);plot(t_lin_conv_vec,lin_conv_vec1,'-s')
lin_conv_vec2=cconv(x_vec,y_vec,2*K-1);plot(t_lin_conv_vec,lin_conv_vec2,'-x')
xlabel('$k$', 'interpreter', 'latex')
ylabel('$\left(x\overline{*}h\right)_{k}$', 'interpreter', 'latex')
title('Linear Convolution');
legend({'FFT','conv','cconv'})

%calculate and plot correlation
subplot(3,2,6)
lin_corr_vec=slow_xcorr(x_vec,y_vec);plot(tau_lin_corr_vec,lin_corr_vec,'-^');hold on;
lin_corr_vec1=ifft(fft(y_vec,2*K-1).*conj(fft(x_vec,2*K-1)));plot(0:2*K-2,lin_corr_vec1,'-o')
lin_corr_vec2=fliplr(conv(x_vec,fliplr(y_vec)));plot(tau_lin_corr_vec,lin_corr_vec2,'-s')
lin_corr_vec3=xcorr(y_vec,x_vec);plot(tau_lin_corr_vec,lin_corr_vec3,'-x');
xlabel('$\kappa$', 'interpreter', 'latex')
ylabel('$\hat{r}_{xy,\kappa}^{\textrm{lin, biased}}$', 'interpreter', 'latex')
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