clearvars
clc
close all

load CantiliverBeamImpactTest    %loads x_rows and y_rows
SNR_db=40;


[P,N]=size(x_rows);
f_s=256;
[T,D_t,D_f]=samplingParameters_fs_N(f_s,N);
t_row=(0:N-1)*D_t;
f_row=(0:N-1)*D_f;

Beta=1.28;  %over sampling facor defined as f_s=2*f_max*Beta
f_max=f_s/2/Beta;
N_f_max=round(N/2/Beta);

x_window_row=ones(1,N);
y_window_row=x_window_row;
x_window_row([1:4,20:end])=0;
y_window_row(1:4)=0;
x_weighted_rows=addNoise(x_rows,db2pow(SNR_db)).*(ones(P,1)*x_window_row);
y_weighted_rows=addNoise(y_rows,db2pow(SNR_db)).*(ones(P,1)*y_window_row);

X_weighted_rows=fft(x_weighted_rows,[],2);
Y_weighted_rows=fft(y_weighted_rows,[],2);

%Plot the first measured x and y signals along with their FFT
figure
subplot(2,2,1)
plot(t_row/T,x_weighted_rows(1,:));
xlim([-0.05,1])
set(gca,'XTickLabel',[]);
ylabel('$x_{1}(t)$','interpreter','latex');

subplot(2,2,2)
semilogy(f_row(1:N_f_max)/f_s,abs(X_weighted_rows(1,1:N_f_max)));
set(gca,'XTickLabel',[]);
ylabel('$\left|X_{1}(f)\right|$','interpreter','latex');

subplot(2,2,3)
plot(t_row/T,y_weighted_rows(1,:));
xlim([-0.05,1])
ylabel('$y_{1}(t)$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');

subplot(2,2,4)
semilogy(f_row(1:N_f_max)/f_s,abs(Y_weighted_rows(1,1:N_f_max)));
ylabel('$\left|Y_{1}(f)\right|$','interpreter','latex');
xlabel('$f/f_{\mathrm{s}}$','interpreter','latex');

%H_raw
H_raw=Y_weighted_rows(1,:)./X_weighted_rows(1,:);

%H_Welch
R_XY=mean(Y_weighted_rows.*conj(X_weighted_rows));
R_XX=mean(abs(X_weighted_rows).^2);
H_Welch=R_XY./R_XX;

%Plot H_raw and H_Welch estimator
figure
ax_mag_h=plot_FRF_mag_phase(f_row(1:N_f_max)/f_s,[H_raw(1:N_f_max);H_Welch(1:N_f_max)],false,[],[],'$f/f_{\mathrm{s}}$');
legend(ax_mag_h,{'$\hat{H}_{\mathrm{raw}}(f)$','$\hat{H}_{\mathrm{Welch}}(f)$'},'interpreter','latex')

export_figure(1,'',{'CantiliverBeamImpactTest_Signals'})
export_figure(2,'',{'CantiliverBeamImpactTest_FRF'})