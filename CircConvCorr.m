clc
close all
clearvars

set(groot,'DefaultAxesColorOrder',[1,0,0])
set(groot,'DefaultAxesLineStyleOrder','-o|-x|-^')

x_vec=[.5,1,1,1,0];
K=length(x_vec);
y_vec=linspace(1,0,K);y_vec=[0,y_vec(1:K-1)];
t_vec=0:K-1;

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
circ_conv_vec=ifft(fft(x_vec).*fft(y_vec));
circ_conv_vec1=cconv(x_vec,y_vec,K);
plot(t_vec,circ_conv_vec,t_vec,circ_conv_vec1);
xlabel('$k$','interpreter','latex')
ylabel('$\left(x * h\right)_{k}$','interpreter','latex')
title('Circular Convolution')
legend(["FFT","cconv"],'Location','southeast')

%calculate and plot correlation
subplot(3,2,6)
circ_corr_vec=ifft(conj(fft(x_vec)).*fft(y_vec));
circ_corr_vec1=fliplr(cconv(x_vec,fliplr(y_vec),K));
circ_corr_vec2=ccorrFunc(x_vec,y_vec);
plot(t_vec,circ_corr_vec,t_vec,circ_corr_vec1,t_vec,circ_corr_vec2);
xlabel('$\kappa$','interpreter','latex')
ylabel('$r_{xy,\kappa}^{\mathrm{circ}}$','interpreter','latex')
title('Circular Correlation')
legend(["FFT","cconv","sum"],'Location','northeast')

export_figure(gcf,'==',"CircConvCorr")

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultAxesLineStyleOrder','remove')
