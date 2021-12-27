clc
close all
clearvars

set(groot,'DefaultAxesLineStyleOrder','o-|x-|^-')
set(groot,'DefaultAxesClipping','off')

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
figure
ax=subplot(4,2,1);
title('$K_{\mathrm{z}}\geq2K-1$','interpreter','latex')
yyaxis left
plot(k_z_vec,x_z_vec)
ylabel("$x_{\mathrm{z},k_{\mathrm{z}}}$",'interpreter','latex')
ylims=ax.YLim;
yyaxis right
plot(k_vec,x_vec)
ylabel("$x_{k}$",'interpreter','latex')
ax.XTick=[0,K-1,K_z-1];
ax.XTickLabel=[];
ax.XAxis.TickLabelInterpreter='latex';
ax.XLim=[0,K_z-1];
[miny,maxy]=bounds([ylims,ax.YLim]);
ylim([miny,maxy]);
ax.XGrid='on';

ax=subplot(4,2,2);
yyaxis left
plot(n_z_fine_vec/N_z_fine,abs(X_z_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1])
hold on
ax.ColorOrderIndex=1;
plot(n_z_vec/N_z,abs(X_z_vec(1:N_z_no_fold)))
hold off
ylabel("$|X_{{\mathrm{z}},n_{\mathrm{z}}}|$",'interpreter','latex')
yyaxis right
plot(n_vec/N,abs(X_vec(1:N_no_fold)))
ylabel("$|X_{n}|$",'interpreter','latex')
ax.XTickLabel=[];
legend(["$|X_{\mathrm{sw}}(f)|$ (DTSTFT)","$|X_{{\mathrm{z}},n_{\mathrm{z}}}|$","$|X_{n}|$"],'interpreter','latex')

ax=subplot(4,2,3);
yyaxis left
plot(k_z_vec,y_z_vec)
ylabel("$y_{\mathrm{z},k_{\mathrm{z}}}$",'interpreter','latex')
ylims=ax.YLim;
yyaxis right
plot(k_vec,y_vec)
ylabel("$y_{k}$",'interpreter','latex')
ax.XTick=[0,K-1,K_z-1];
ax.XTickLabel=[];
ax.XAxis.TickLabelInterpreter='latex';
ax.XLim=[0,K_z-1];
[miny,maxy]=bounds([ylims,ax.YLim]);
ylim([miny,maxy]);
ax.XGrid='on';

ax=subplot(4,2,4);
yyaxis left
plot(n_z_fine_vec/N_z_fine,abs(Y_z_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1])
hold on
ax.ColorOrderIndex=1;
plot(n_z_vec/N_z,abs(Y_z_vec(1:N_z_no_fold)))
hold off
ylabel("$|Y_{{\mathrm{z}},n_{\mathrm{z}}}|$",'interpreter','latex')
yyaxis right
plot(n_vec/N,abs(Y_vec(1:N_no_fold)))
ylabel("$|Y_{n}|$",'interpreter','latex')
ax.XTickLabel=[];
legend(["$|Y_{\mathrm{sw}}(f)|$ (DTSTFT)","$|Y_{{\mathrm{z}},n_{\mathrm{z}}}|$","$|Y_{n}|$"],'interpreter','latex')

%Convolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
circ_Conv_vec=X_vec.*Y_vec*D_t;
lin_Conv_vec=X_z_vec.*Y_z_vec*D_t;
DTST_Conv_fine_vec=X_z_fine_vec.*Y_z_fine_vec*D_t;

circ_conv_vec=cconv(x_vec,y_vec,K);
circ_conv_vec1=ifft(circ_Conv_vec);
lin_conv_vec=conv(x_vec,y_vec)*D_t;
lin_conv_vec1=cconv(x_vec,y_vec,N_z)*D_t;
lin_conv_vec2=ifft(lin_Conv_vec);

ax=subplot(4,2,5);
yyaxis left
plot(kappa_lin_conv_vec,[lin_conv_vec;lin_conv_vec1;lin_conv_vec2])
ylabel('$c_{k}^{\mathrm{lin}} \equiv \left(x\overline{*}y\right)_{k}=\left(x_{\mathrm{z}}\bigotimes y_{\mathrm{z}}\right)_{k_{\mathrm{z}}=k}$','interpreter','latex')
ylims=ax.YLim;
yyaxis right
plot(k_vec,[circ_conv_vec;circ_conv_vec1]);
ylabel('$c_{k}^{\mathrm{circ}} \equiv \left(x \bigotimes y\right)_{k}$','interpreter','latex')
ax.XTick=[0,K-1,K_z-1];
ax.XTickLabel=["$0$","$K-1$","$K_{\mathrm{z}}-1=2K-2$"];
ax.XAxis.TickLabelInterpreter='latex';
ax.XLim=[0,K_z-1];
legend(["$c_{k}^{\mathrm{lin}}$ (conv)","$\hphantom{c_{k}^{\mathrm{lin}}}$ (cconv)","$\hphantom{c_{k}^{\mathrm{lin}}}$ (FFT)","$c_{k}^{\mathrm{circ}}$ (cconv)","$\hphantom{c_{k}^{\mathrm{circ}}}$ (FFT)"],'interpreter','latex')
ax.XGrid='on';
[miny,maxy]=bounds([ylims,ax.YLim]);
ylim([miny,maxy]);
ax.XGrid='on';

ax=subplot(4,2,6);
yyaxis left
plot(n_z_fine_vec/N_z_fine,abs(DTST_Conv_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1])
hold on
ax.ColorOrderIndex=1;
plot(n_z_vec/N_z,abs(lin_Conv_vec(1:N_z_no_fold)))
hold off
ylabel("$C_{n}^{\mathrm{lin}}=\left|X_{\mathrm{z},n_{\mathrm{z}}}\,Y_{\mathrm{z},n_{\mathrm{z}}}\right|_{n_{z}=n}$",'interpreter','latex')
yyaxis right
plot(n_vec/N,abs(circ_Conv_vec(1:N_no_fold)))
ylabel("$C_{n}^{\mathrm{circ}} = |X_{n} . Y_{n}|$",'interpreter','latex')
legend(["$|C_{\mathrm{sw}}(f)|$ (DTSTFT)","$|C_{n}^{\mathrm{lin}}|$","$|C_{n}^{\mathrm{circ}}|$"],'interpreter','latex')
ax.XTickLabel=[];

%calculate and plot correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
circ_Corr_vec=conj(X_vec).*Y_vec;
lin_Corr_vec=Y_z_vec.*conj(X_z_vec)*D_t;
lin_Corr_fine_vec=Y_z_fine_vec.*conj(X_z_fine_vec)*D_t;

circ_corr_vec=ccorrFunc(x_vec,y_vec);
circ_corr_vec1=fliplr(cconv(x_vec,fliplr(y_vec),K));
circ_corr_vec2=ifft(circ_Corr_vec);

lin_corr_vec=slow_xcorr(x_vec,y_vec)*D_t;
lin_corr_vec1=fliplr(conv(x_vec,fliplr(y_vec)))*D_t;
lin_corr_vec2=xcorr(y_vec,x_vec)*D_t;
lin_corr_vec3=fftshift(ifft(lin_Corr_vec));

ax=subplot(4,2,7);
yyaxis left
plot(kappa_lin_corr_vec,[lin_corr_vec;lin_corr_vec1;lin_corr_vec2;lin_corr_vec3]);
ylims=ax.YLim;
ylabel('$r_{xy,\kappa}^{\mathrm{lin}}=r_{x_{\mathrm{z}}y_{\mathrm{z}},\kappa_{\mathrm{z}}=\kappa}^{\mathrm{circ}}$','interpreter','latex')
yyaxis right
plot(k_vec,[circ_corr_vec;circ_corr_vec1;circ_corr_vec2]);
ylabel('$r_{xy,\kappa}^{\mathrm{circ}}$','interpreter','latex')
ax.XTick=[-(K-1),0,K-1];
ax.XTickLabel=["$-(K-1)$","$0$","$K-1$"];
ax.XAxis.TickLabelInterpreter='latex';
ax.XLim=[-(K-1),K-1];
xlabel('$k,\quad k_{\mathrm{z}},\quad \kappa$','interpreter','latex')
legend(["$r_{xy,\kappa}^{\mathrm{lin}}$ $(\sum)$","$\hphantom{\ensuremath{r_{xy,\kappa}^{\mathrm{lin}}}}$ (cconv)","$\hphantom{\ensuremath{r_{xy,\kappa}^{\mathrm{lin}}}}$ (xcorr)","$\hphantom{\ensuremath{r_{xy,\kappa}^{\mathrm{lin}}}}$ (FFT)","$r_{xy,\kappa}^{\mathrm{circ}}$ ($\sum$)","$\hphantom{r_{xy,\kappa}^{\mathrm{circ}}}$ (cconv)","$\hphantom{r_{xy,\kappa}^{\mathrm{circ}}}$ (FFT)"],'interpreter','latex','Location','northwest')
ax.XGrid='on';
[miny,maxy]=bounds([ylims,ax.YLim]);
ylim([miny,maxy]);
ax.XGrid='on';

ax=subplot(4,2,8);
yyaxis left
plot(n_z_fine_vec/N_z_fine,abs(lin_Corr_fine_vec(1:N_z_fine_no_fold)),'-','Color',.4*[1,1,1])
hold on
ax.ColorOrderIndex=1;
plot(n_z_vec/N_z,abs(lin_Corr_vec(1:N_z_no_fold)));
ylabel("$R_{XY,n}^{\mathrm{lin}}=R_{X_{\mathrm{z}}Y_{\mathrm{z}},n_{\mathrm{z}}=n}^{\mathrm{circ}}$",'interpreter','latex')
yyaxis right
plot(n_vec/N,abs(circ_Corr_vec(1:N_no_fold)));
ylabel("$R_{XY,n}^{\mathrm{circ}} \equiv Y_{n}\,X_{n}^{*}$",'interpreter','latex')
legend(["$|R_{XY,\mathrm{sw}}(f)|$ (DTSTFT)","$|R_{XY,n}^{\mathrm{lin}}|$","$|R_{XY,n}^{\mathrm{circ}}|$"],'interpreter','latex')
xlabel('$n/N,\quad n_{\mathrm{z}}/N_{\mathrm{z}},\quad f/f_{\mathrm{s}}$','interpreter','latex');

%Additional optimization of the axes for correct comparison with the correlation curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos=get(subplot(4,2,1),'Position');
width=pos(3);
for ii=[1,3,5]
    pos=get(subplot(4,2,ii),'Position');
    pos(1)=pos(1)+width/3;
    pos(3)=width*2/3;
    set(subplot(4,2,ii),'Position',pos)
end

pos=get(subplot(4,2,7),'Position');
pos(3)=width*2/3;
set(subplot(4,2,7),'Position',pos)

for ii=2:2:8
    pos=get(subplot(4,2,ii),'Position');
    pos(3)=width*2/3;
    set(subplot(4,2,ii),'Position',pos)
end

export_figure(gcf,'==',"DigitalConvCorr")

set(groot,'DefaultAxesLineStyleOrder','remove')
set(groot,'DefaultAxesClipping','remove')