clc
% close all
clearvars

set(groot,'DefaultAxesColorOrder',[1,0,0])
set(groot,'DefaultAxesLineStyleOrder','-o|-x|-^|-s')

x_vec=[1,.75,.5,.25,.25];
K=length(x_vec);
y_vec=hamming(5,"periodic").';
D_t=1;
k_vec=0:K-1;

K_z=2*K-1;
kappa_lin_conv_vec=0:K_z-1;
kappa_lin_corr_vec=-(K-1):K-1;

N=K;
N_z=K_z;
[~,~,~,N_z_no_fold]=samplingParameters_D_t_N(D_t,N_z);
X_z_vec=fft(x_vec,N_z);
Y_z_vec=fft(y_vec,N_z);
n_z_vec=0:N_z_no_fold-1;

N_z_fine=20*N;
[~,~,~,N_z_fine_no_fold]=samplingParameters_D_t_N(D_t,N_z_fine);
X_z_fine_vec=fft(x_vec,N_z_fine);
Y_z_fine_vec=fft(y_vec,N_z_fine);
n_z_fine_vec=0:N_z_fine_no_fold-1;

%plot x and y signals
figures=[figure,figure];
titles="Linear "+["Convolution","Correlation"];
for ii=1:2
    figure(figures(ii))

    subplot(3,2,1);
    plot(k_vec/K_z,x_vec)
    ylabel('$x_{\mathrm{z},k_{\mathrm{z}}=k}=x_{k}\;\forall\:0\leq k_{\mathrm{z}}\leq K-1$','interpreter','latex')
    xlabel('$k_{\mathrm{z}}/K_{\mathrm{z}}$','interpreter','latex');
    
    subplot(3,2,2);
    plot(n_z_vec/N_z,abs(X_z_vec(1:N_z_no_fold)))
    hold on,plot(n_z_fine_vec/N_z_fine,abs(X_z_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1],'DisplayName','$|X_{\mathrm{sw}}(f)|$ (DTSTFT)');hold off
    ylabel('$|X_{{\mathrm{z}},n_{\mathrm{z}}}|$','interpreter','latex')
    xlabel('$n_{\mathrm{z}}/N_{\mathrm{z}}$','interpreter','latex');

    subplot(3,2,3);
    plot(k_vec/K_z,y_vec)
    ylabel('$y_{\mathrm{z},k_{\mathrm{z}}=k}=y_{k}\;\forall\:0\leq k_{\mathrm{z}}\leq K-1$','interpreter','latex')
    xlabel('$k_{\mathrm{z}}/K_{\mathrm{z}}$','interpreter','latex');

    subplot(3,2,4);
    plot(n_z_vec/N_z,abs(Y_z_vec(1:N_z_no_fold)))
    hold on,plot(n_z_fine_vec/N_z_fine,abs(Y_z_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1],'DisplayName','$|Y_{\mathrm{sw}}(f)|$ (DTSTFT)');hold off
    ylabel('$|Y_{{\mathrm{z}},n_{\mathrm{z}}}|$','interpreter','latex')
    xlabel('$n_{\mathrm{z}}/N_{\mathrm{z}}$','interpreter','latex');

    sgtitle(titles(ii))
end

%calculate and plot convolution
lin_Conv_vec=X_z_vec.*Y_z_vec*D_t;
lin_Conv_fine_vec=X_z_fine_vec.*Y_z_fine_vec*D_t;

figure(figures(1))
subplot(3,2,5)
lin_conv_vec=conv(x_vec,y_vec)*D_t;
lin_conv_vec1=cconv(x_vec,y_vec,N_z)*D_t;
lin_conv_vec2=ifft(lin_Conv_vec);
plot(kappa_lin_conv_vec/K_z,[lin_conv_vec;lin_conv_vec1;lin_conv_vec2])
ylabel('$\left(x\overline{*}h\right)_{k}=\left(x_{\mathrm{z}}\bigotimes h_{\mathrm{z}}\right)_{k_{\mathrm{z}}=k} \; \forall \: 0\leq k\leq2K-2$','interpreter','latex')
xlabel('$k_{\mathrm{z}}/K_{\mathrm{z}}$','interpreter','latex');
legend(["conv","cconv","FFT"])

subplot(3,2,6)
plot(n_z_vec/N_z,abs(lin_Conv_vec(1:N_z_no_fold)))
hold on,plot(n_z_fine_vec/N_z_fine,abs(lin_Conv_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1],'DisplayName','$|Y_{\mathrm{sw}}(f)|$ (DTSTFT)');hold off
ylabel('$\left|X_{\mathrm{z},n_{\mathrm{z}}}\,Y_{\mathrm{z},n_{\mathrm{z}}}\right|$','interpreter','latex')
xlabel('$n_{\mathrm{z}}/N_{\mathrm{z}}$','interpreter','latex');

%calculate and plot correlation
lin_Corr_vec=Y_z_vec.*conj(X_z_vec)*D_t;
lin_Corr_fine_vec=Y_z_fine_vec.*conj(X_z_fine_vec)*D_t;

figure(figures(2))
subplot(3,2,5)
lin_corr_vec=slow_xcorr(x_vec,y_vec)*D_t;
lin_corr_vec1=fliplr(conv(x_vec,fliplr(y_vec)))*D_t;
lin_corr_vec2=xcorr(y_vec,x_vec)*D_t;
lin_corr_vec3=fftshift(ifft(lin_Corr_vec));
plot(kappa_lin_corr_vec/K_z,[lin_corr_vec;lin_corr_vec1;lin_corr_vec2;lin_corr_vec3]);
ylabel('$r_{xy,\kappa}^{\mathrm{lin}}=r_{x_{\mathrm{z}}y_{\mathrm{z}},\kappa_{\mathrm{z}}=\kappa}^{\mathrm{circ}} \; \forall \: -(K-1)\leq\kappa\leq K-1$','interpreter','latex')
xlabel('$\kappa /K_{\mathrm{z}}$','interpreter','latex')
legend(["$\sum$","cconv","xcorr","FFT"],'interpreter','latex','Location','northwest')

subplot(3,2,6)
plot(n_z_vec/N_z,abs(lin_Corr_vec(1:N_z_no_fold)));
hold on,plot(n_z_fine_vec/N_z_fine,abs(lin_Corr_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1],'DisplayName','$|Y_{\mathrm{sw}}(f)|$ (DTSTFT)');hold off
ylabel('$Y_{\mathrm{z},n_{\mathrm{z}}}\,\left(X_{\mathrm{z},n_{\mathrm{z}}}\right)^{*}\equiv R_{X_{\mathrm{z}}Y_{\mathrm{z}},n_{\mathrm{z}}}^{\mathrm{circ}}$','interpreter','latex')
xlabel('$n_{\mathrm{z}}$','interpreter','latex')

%Additional optimization of the axes for correct comparison with the correlation curves
figure(figures(1))
for ii=[1,3]
    pos=get(subplot(3,2,ii),'Position');
    pos(3)=pos(3)/2;
    set(subplot(3,2,ii),'Position',pos)
    xlim([-inf,0.5])
end

figure(figures(2))
for ii=[1,3]
    pos=get(subplot(3,2,ii),'Position');
    pos(1)=pos(1)+pos(3)/2;
    pos(3)=pos(3)/2;
    set(subplot(3,2,ii),'Position',pos)
    xlim([-inf,0.5])
end

export_figure(gcf,'==',"LinearConvCorr")

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultAxesLineStyleOrder','remove')