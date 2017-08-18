function compare_H_estimators()
clc
close all
set(groot,'DefaultLineLineWidth',1);

load CantiliverBeamImpactTest    %loads x_rows and y_rows

[Q,N]=size(x_rows);
f_s=256;
[T,d_t,d_f]=samplingParameters_fs_N(f_s,N);

Beta=1.28;  %over sampling facor defined as f_s=2*f_max*Beta
f_max=f_s/2/Beta;
N_f_max=round(N/2/Beta);

f=[0:N-1]*d_f;
t=[0:N-1]*d_t;

X_rows=fft(x_rows,[],2);
Y_rows=fft(y_rows,[],2);

%calculate H_raw
H_raw=Y_rows(1,:)./X_rows(1,:);
h_raw=ifft(H_raw);

%H1 estimator
R_XY=mean(Y_rows.*conj(X_rows));
R_XX=mean(X_rows.*conj(X_rows));
H_1=R_XY./R_XX;

%Plot the first measured x and y along with their FFT
figure
subplot(3,2,1)
plot(t, x_rows(1,:));
xlim([-0.5,T])
ylabel('$x_{1}(t)$', 'interpreter', 'latex');

subplot(3,2,2)
semilogy(f(1:N_f_max), abs(X_rows(1,1:N_f_max)));
axis tight
ylabel('$\left|X_{1}(f)\right|$', 'interpreter', 'latex');

subplot(3,2,3)
plot(t, y_rows(1,:));
xlim([-0.5,T])
ylabel('$y_{1}(t)$', 'interpreter', 'latex');

subplot(3,2,4)
semilogy(f(1:N_f_max), abs(Y_rows(1,1:N_f_max)));
axis tight
ylabel('$\left|Y_{1}(f)\right|$', 'interpreter', 'latex');

%Plot h_raw & H_raw
subplot(3,2,5);
plot(t, h_raw);
xlim([-0.5,T])
xlabel('$t$ (s)', 'interpreter', 'latex');
ylabel('$h_{\mathrm{raw}}(t)$', 'interpreter', 'latex');

subplot(3,2,6);
semilogy(f(1:N_f_max),abs(H_raw(1:N_f_max)));
axis tight
xlabel('$f$ (Hz)', 'interpreter', 'latex');
ylabel('$\left|H_{\mathrm{raw}}(f)\right|$', 'interpreter', 'latex');


%Plot H_raw beside H_1 estimator
figure
ax_mag_h=subplot(4,1,[1:3]);hold on;axis tight
ax_phase_h=subplot(4,1,4);hold on;axis tight
plot_FRF_mag_phase(f(1:N_f_max),H_raw(1:N_f_max),false,ax_mag_h,ax_phase_h,[]);
plot_FRF_mag_phase(f(1:N_f_max),H_1(1:N_f_max),false,ax_mag_h,ax_phase_h,[]);
legend(ax_mag_h,{'$\left|H_{\mathrm{raw}}(f)\right|$','$\left|\hat{H}_{1}(f)\right|$'}, 'interpreter', 'latex')
legend(ax_phase_h,{'$\angle \left(H_{\mathrm{raw}}(f)\right)$ (rad)','$\angle \left(\hat{H}_{1}(f)\right)$ (rad)'}, 'interpreter', 'latex')



set(groot,'DefaultLineLineWidth','remove');

export_figure(1,'||',{'CantiliverBeamImpactTest_Signals'})
export_figure(2,'',{'CantiliverBeamImpactTest_FRF'})