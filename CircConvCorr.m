clc
close all
clearvars

set(groot,'DefaultAxesColorOrder',[1,0,0])
set(groot,'DefaultAxesLineStyleOrder','-o|-x|-^')

x_vec=[1,.75,.5,.25,.25];
K=length(x_vec);
y_vec=hamming(5,'periodic').';
D_t=1;
k_vec=0:K-1;

N=K;
[~,~,~,N_no_fold]=samplingParameters_D_t_N(D_t,N);
X_vec=fft(x_vec);
Y_vec=fft(y_vec);
n_vec=0:N_no_fold-1;

N_z_fine=20*N;
[~,~,~,N_z_fine_no_fold]=samplingParameters_D_t_N(D_t,N_z_fine);
X_z_fine_vec=fft(x_vec,N_z_fine);
Y_z_fine_vec=fft(y_vec,N_z_fine);
n_z_fine_vec=0:N_z_fine_no_fold-1;

%plot x and y signals
figures=[figure,figure];
titles="Circular "+["Convolution","Correlation"];
for ii=1:2
    figure(figures(ii))

    subplot(3,2,1);
    plot(k_vec/K,x_vec)
    ylabel('$x_{k}$','interpreter','latex')
    xlabel('$k/K$','interpreter','latex');
    xlim([0,1])
    
    subplot(3,2,2);
    plot(n_vec/N,abs(X_vec(1:N_no_fold)))
    hold on,plot(n_z_fine_vec/N_z_fine,abs(X_z_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1],'DisplayName','$|X_{\mathrm{sw}}(f)|$ (DTSTFT)');hold off
    ylabel('$|X_{n}|$','interpreter','latex')
    xlabel('$n/N$','interpreter','latex');
    xlim([0,.5])

    subplot(3,2,3);
    plot(k_vec/K,y_vec)
    ylabel('$y_{k}$','interpreter','latex')
    xlabel('$k/K$','interpreter','latex');
    xlim([0,1])
    
    subplot(3,2,4);
    plot(n_vec/N,abs(Y_vec(1:N_no_fold)))
    hold on,plot(n_z_fine_vec/N_z_fine,abs(Y_z_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1],'DisplayName','$|Y_{\mathrm{sw}}(f)|$ (DTSTFT)');hold off
    ylabel('$|Y_{n}|$','interpreter','latex')
    xlabel('$n/N$','interpreter','latex');
    xlim([0,.5])

    sgtitle(titles(ii))
end

%calculate and plot convolution
circ_Conv_vec=X_vec.*Y_vec*D_t;
circ_Conv_fine_vec=X_z_fine_vec.*Y_z_fine_vec*D_t;

figure(figures(1))
subplot(3,2,5)
circ_conv_vec=ifft(circ_Conv_vec);
circ_conv_vec1=cconv(x_vec,y_vec,K);
plot(k_vec/K,[circ_conv_vec;circ_conv_vec1]);
xlim([0,1])
xlabel('$k/K$','interpreter','latex')
ylabel('$\left(x \bigotimes h\right)_{k}$','interpreter','latex')
legend(["FFT","cconv"],'Location','southeast')

subplot(3,2,6)
plot(n_vec/N,abs(circ_Conv_vec(1:N_no_fold)))
hold on,plot(n_z_fine_vec/N_z_fine,abs(circ_Conv_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1],'DisplayName','$|Y_{\mathrm{sw}}(f)|$ (DTSTFT)');hold off
xlim([0,.5])
xlabel('$n/N$','interpreter','latex')
ylabel('$|X_{n} . H_{n}|$','interpreter','latex')

%calculate and plot correlation
circ_Corr_vec=conj(X_vec).*Y_vec;
circ_Corr_fine_vec=conj(X_z_fine_vec).*Y_z_fine_vec;

figure(figures(2))
subplot(3,2,5)
circ_corr_vec=ifft(circ_Corr_vec);
circ_corr_vec1=fliplr(cconv(x_vec,fliplr(y_vec),K));
circ_corr_vec2=ccorrFunc(x_vec,y_vec);
plot(k_vec/K,[circ_corr_vec;circ_corr_vec1;circ_corr_vec2]);
xlim([0,1])
ylabel('$r_{xy,\kappa}^{\mathrm{circ}}$','interpreter','latex')
xlabel('$\kappa /K$','interpreter','latex')
legend(["FFT","cconv","$\sum$"],'interpreter','latex','Location','northeast')

subplot(3,2,6)
plot(n_vec/N,abs(circ_Corr_vec(1:N_no_fold)));
hold on,plot(n_z_fine_vec/N_z_fine,abs(circ_Corr_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1],'DisplayName','$|Y_{\mathrm{sw}}(f)|$ (DTSTFT)');hold off
xlim([0,.5])
ylabel('$R_{XY,n}^{\mathrm{circ}}$','interpreter','latex')
xlabel('$n/N$','interpreter','latex')

export_figure(gcf,'==',"CircConvCorr")

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultAxesLineStyleOrder','remove')
