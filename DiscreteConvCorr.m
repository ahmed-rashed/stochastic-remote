clc
% close all
clearvars

% set(groot,'DefaultAxesColorOrder',[1,0,0])
set(groot,'DefaultAxesLineStyleOrder','-o|-x|-^|-s')

x_vec=[1,.75,.5,.25,.25];
K=length(x_vec);
y_vec=hamming(5,"periodic").';

D_t=1;
k_vec=0:K-1;
N=K;
[~,~,~,N_no_fold]=samplingParameters_D_t_N(D_t,N);
X_vec=fft(x_vec);
Y_vec=fft(y_vec);
n_vec=0:N_no_fold-1;

K_z=2*K-1;
k_z_vec=0:K_z-1;
x_z_vec=[x_vec,zeros(1,K_z-K)];
y_z_vec=[y_vec,zeros(1,K_z-K)];
kappa_lin_conv_vec=0:K_z-1;
kappa_lin_corr_vec=-(K-1):K-1;
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
titles=["Convolution","Correlation"];
for ii=1:2
    figure(figures(ii))

    ax=subplot(3,2,1);
    title('$K_{\mathrm{z}}\geq2K-1$','interpreter','latex')
    hold on
    plot(k_z_vec,x_z_vec)
    plot(k_vec,x_vec)
    hold off
    legend(["$x_{\mathrm{z},k_{\mathrm{z}}}$","$x_{k}$"],'interpreter','latex')
    ax.XTick=[0,K-1,K_z-1];
    ax.XTickLabel=["$0$","$K-1$","$K_{\mathrm{z}}-1$"];
    ax.XAxis.TickLabelInterpreter='latex';
    ax.XLim=[0,K_z-1];
    ax.XGrid='on';
    
    ax=subplot(3,2,2);
    hold on
    plot(n_z_fine_vec/N_z_fine,abs(X_z_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1])
    plot(n_z_vec/N_z,abs(X_z_vec(1:N_z_no_fold)))
    plot(n_vec/N,abs(X_vec(1:N_no_fold)))
    hold off
    ax.XTickLabel=[];
    legend(["$|X_{\mathrm{sw}}(f)|$ (DTSTFT)","$|X_{{\mathrm{z}},n_{\mathrm{z}}}|$ (DTSTFT)","$|X_{n}|$ (DFT)"],'interpreter','latex')

    ax=subplot(3,2,3);
    hold on
    plot(k_z_vec,y_z_vec)
    plot(k_vec,y_vec)
    hold off
    legend(["$y_{\mathrm{z},k_{\mathrm{z}}}$","$y_{k}$"],'interpreter','latex')
    ax.XTick=[0,K-1,K_z-1];
    ax.XTickLabel=["$0$","$K-1$","$K_{\mathrm{z}}-1$"];
    ax.XAxis.TickLabelInterpreter='latex';
    ax.XLim=[0,K_z-1];
    ax.XGrid='on';

    ax=subplot(3,2,4);
    hold on
    plot(n_z_fine_vec/N_z_fine,abs(Y_z_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1])
    plot(n_z_vec/N_z,abs(Y_z_vec(1:N_z_no_fold)))
    plot(n_vec/N,abs(Y_vec(1:N_no_fold)))
    hold off
    ax.XTickLabel=[];
    legend(["$|Y_{\mathrm{sw}}(f)|$ (DTSTFT)","$|Y_{{\mathrm{z}},n_{\mathrm{z}}}|$ (DTSTFT)","$|Y_{n}|$ (DFT)"],'interpreter','latex')

    sgtitle(titles(ii))
end

%calculate and plot convolution
circ_Conv_vec=X_vec.*Y_vec*D_t;
lin_Conv_vec=X_z_vec.*Y_z_vec*D_t;
DTST_Conv_fine_vec=X_z_fine_vec.*Y_z_fine_vec*D_t;

figure(figures(1))
circ_conv_vec=ifft(circ_Conv_vec);
circ_conv_vec1=cconv(x_vec,y_vec,K);
lin_conv_vec=conv(x_vec,y_vec)*D_t;
lin_conv_vec1=cconv(x_vec,y_vec,N_z)*D_t;
lin_conv_vec2=ifft(lin_Conv_vec);
ax=subplot(3,2,5);
hold on
plot(k_vec,[circ_conv_vec;circ_conv_vec1]);
plot(kappa_lin_conv_vec,[lin_conv_vec;lin_conv_vec1;lin_conv_vec2])
hold off
ax.XTick=[0,K-1,K_z-1];
ax.XTickLabel=["$0$","$K-1$","$K_{\mathrm{z}}-1=2K-2$"];
ax.XAxis.TickLabelInterpreter='latex';
ax.XLim=[0,K_z-1];
ylabel('$\left(x\overline{*}h\right)_{k}=\left(x_{\mathrm{z}}\bigotimes h_{\mathrm{z}}\right)_{k_{\mathrm{z}}=k}$','interpreter','latex')
xlabel('$k,\quad k_{\mathrm{z}}$','interpreter','latex');
legend(["FFT","cconv","conv","cconv","FFT"])
ax.XGrid='on';

subplot(3,2,6)
hold on
plot(n_vec/N,abs(circ_Conv_vec(1:N_no_fold)))
plot(n_z_vec/N_z,abs(lin_Conv_vec(1:N_z_no_fold)))
plot(n_z_fine_vec/N_z_fine,abs(DTST_Conv_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1])
hold off
ylabel('$\left|X_{\mathrm{z},n_{\mathrm{z}}}\,Y_{\mathrm{z},n_{\mathrm{z}}}\right|$','interpreter','latex')
xlabel('$n/N,\quad n_{\mathrm{z}}/N_{\mathrm{z}},\quad f/f_{\mathrm{s}}$','interpreter','latex');

%calculate and plot correlation
circ_Corr_vec=conj(X_vec).*Y_vec;
lin_Corr_vec=Y_z_vec.*conj(X_z_vec)*D_t;
lin_Corr_fine_vec=Y_z_fine_vec.*conj(X_z_fine_vec)*D_t;

figure(figures(2))
circ_corr_vec=ifft(circ_Corr_vec);
circ_corr_vec1=fliplr(cconv(x_vec,fliplr(y_vec),K));
circ_corr_vec2=ccorrFunc(x_vec,y_vec);

lin_corr_vec=slow_xcorr(x_vec,y_vec)*D_t;
lin_corr_vec1=fliplr(conv(x_vec,fliplr(y_vec)))*D_t;
lin_corr_vec2=xcorr(y_vec,x_vec)*D_t;
lin_corr_vec3=fftshift(ifft(lin_Corr_vec));
ax=subplot(3,2,5);
hold on
plot(k_vec,[circ_corr_vec;circ_corr_vec1;circ_corr_vec2]);
plot(kappa_lin_corr_vec,[lin_corr_vec;lin_corr_vec1;lin_corr_vec2;lin_corr_vec3]);
hold off
ax.XTick=[-(K-1),0,K-1];
ax.XTickLabel=["$-(K-1)$","$0$","$K-1$"];
ax.XAxis.TickLabelInterpreter='latex';
ax.XLim=[-(K-1),K-1];
ylabel('$r_{xy,\kappa}^{\mathrm{lin}}=r_{x_{\mathrm{z}}y_{\mathrm{z}},\kappa_{\mathrm{z}}=\kappa}^{\mathrm{circ}}$','interpreter','latex')
xlabel('$k,\quad k_{\mathrm{z}},\quad \kappa$','interpreter','latex')
legend(["FFT","cconv","$\sum$","$\sum$","cconv","xcorr","FFT"],'interpreter','latex','Location','northwest')
ax.XGrid='on';

subplot(3,2,6)
hold on
plot(n_vec/N,abs(circ_Corr_vec(1:N_no_fold)));
plot(n_z_vec/N_z,abs(lin_Corr_vec(1:N_z_no_fold)));
plot(n_z_fine_vec/N_z_fine,abs(lin_Corr_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1])
hold off
ylabel('$Y_{\mathrm{z},n_{\mathrm{z}}}\,\left(X_{\mathrm{z},n_{\mathrm{z}}}\right)^{*}\equiv R_{X_{\mathrm{z}}Y_{\mathrm{z}},n_{\mathrm{z}}}^{\mathrm{circ}}$','interpreter','latex')
xlabel('$n/N,\quad n_{\mathrm{z}}/N_{\mathrm{z}},\quad f/f_{\mathrm{s}}$','interpreter','latex');

%Additional optimization of the axes for correct comparison with the correlation curves
% figure(figures(1))
% for ii=[1,3]
%     pos=get(subplot(3,2,ii),'Position');
%     pos(3)=pos(3)/2;
%     set(subplot(3,2,ii),'Position',pos)
% end

figure(figures(2))
for ii=[1,3]
    pos=get(subplot(3,2,ii),'Position');
    pos(1)=pos(1)+pos(3)/3;
    pos(3)=pos(3)*2/3;
    set(subplot(3,2,ii),'Position',pos)
end

pos=get(subplot(3,2,5),'Position');
pos(3)=pos(3)*2/3;
set(subplot(3,2,5),'Position',pos)

% export_figure(gcf,'==',"LinearConvCorr")

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultAxesLineStyleOrder','remove')