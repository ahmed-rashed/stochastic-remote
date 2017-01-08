function compare_H_estimators_Comprehensive()
clc
close all
set(groot,'DefaultLineLineWidth',1);

load CantiliverBeamImpactTest    %loads x_rows and y_rows

N_avg=size(x_rows,1);
N=size(x_rows,2);
f_s=256;
[T,d_t,d_f]=samplingParameters_fs_N(f_s,N);

Beta=1.28;  %over sampling facor defined as f_s=2*f_max*Beta
f_max=f_s/2/Beta;
N_f_max=round(N/2/Beta);

f=[0:N-1]*d_f;
t=[0:N-1]*d_t;

X_rows=fft(x_rows,[],2);
Y_rows=fft(y_rows,[],2);

x_fig=figure;
y_fig=figure;
for n=1:N_avg
    figure(x_fig);h_x=subplot(N_avg,2,2*n-1);
    plot(t, x_rows(n,:));
    axis tight
    xlim([-0.5,T])
    ylabel(['$x_{',int2str(n),'}(t)$'], 'interpreter', 'latex');
    
    figure(x_fig);h_y=subplot(N_avg,2,2*n);
    semilogy(f(1:N_f_max), abs(X_rows(n,1:N_f_max)));
    axis tight
    ylabel(['$\left|X_{',int2str(n),'}(f)\right|$'], 'interpreter', 'latex');

    if n~=N_avg
        set(h_x,'XTickLabel',[]);
        set(h_y,'XTickLabel',[]);
    else
        xlabel(h_x,'$t$ (s)', 'interpreter', 'latex');
        xlabel(h_y,'$f$ (Hz)', 'interpreter', 'latex');
    end
    
    figure(y_fig);h_x=subplot(N_avg,2,2*n-1);
    plot(t, y_rows(n,:));
    axis tight
    xlim([-0.5,T])
    ylabel(['$y_{',int2str(n),'}(t)$'], 'interpreter', 'latex');
    
    figure(y_fig);h_y=subplot(N_avg,2,2*n);
    semilogy(f(1:N_f_max), abs(Y_rows(n,1:N_f_max)));
    axis tight
    ylabel(['$\left|Y_{',int2str(n),'}(f)\right|$'], 'interpreter', 'latex');

    if n~=N_avg
        set(h_x,'XTickLabel',[]);
        set(h_y,'XTickLabel',[]);
    else
        xlabel(h_x,'$t$ (s)', 'interpreter', 'latex');
        xlabel(h_y,'$f$ (Hz)', 'interpreter', 'latex');
    end
end

figure
subplot(3,2,1)
plot(t, x_rows(1,:));
xlabel('$t$ (s)', 'interpreter', 'latex');
xlim([-0.5,T])
ylabel('$x_{n}(t)$', 'interpreter', 'latex');

subplot(3,2,2)
semilogy(f(1:N_f_max), abs(X_rows(1,1:N_f_max)));
xlabel('$f$ (Hz)', 'interpreter', 'latex');
ylabel('$\left|X_{n}(f)\right|$', 'interpreter', 'latex');

subplot(3,2,3)
plot(t, y_rows(1,:));
xlim([-0.5,T])
xlabel('$t$ (s)', 'interpreter', 'latex');
ylabel('$y_{n}(t)$', 'interpreter', 'latex');

subplot(3,2,4)
semilogy(f(1:N_f_max), abs(Y_rows(1,1:N_f_max)));
xlabel('$f$ (Hz)', 'interpreter', 'latex');
ylabel('$\left|Y_{n}(f)\right|$', 'interpreter', 'latex');

%Raw H estimator
H_raw=Y_rows(1,:)./X_rows(1,:);
h_raw=ifft(H_raw);

subplot(3,2,5);
plot(t, h_raw);
xlim([-0.5,T])
xlabel('$t$ (s)', 'interpreter', 'latex');
ylabel('$h_{\mathrm{raw}}(t)$', 'interpreter', 'latex');

subplot(3,2,6);
semilogy(f(1:N_f_max),abs(H_raw(1:N_f_max)));
xlabel('$f$ (Hz)', 'interpreter', 'latex');
ylabel('$H_{\mathrm{raw}}$', 'interpreter', 'latex');


%H1 estimator
R_XY=mean(Y_rows.*conj(X_rows));
R_XX=mean(X_rows.*conj(X_rows));
R_YY=mean(Y_rows.*conj(Y_rows));

H_1=R_XY./R_XX;
H_2=R_YY./conj(R_XY);

figure
ax_mag_h=subplot(4,3,[1:3:7]);
ax_phase_h=subplot(4,3,10);
plot_FRF_mag_phase(f(1:N_f_max),H_raw(1:N_f_max),false,ax_mag_h,ax_phase_h,[],'H_{\mathrm{raw}}');

% ax_mag_h=subplot(4,3,[1:3:7]+1);
% ax_phase_h=subplot(4,3,10+1);
hold(ax_mag_h,'all')
hold(ax_phase_h,'all')
plot_FRF_mag_phase(f(1:N_f_max),H_1(1:N_f_max),false,ax_mag_h,ax_phase_h,[],'\hat{H}_{1}');

% ax_mag_h=subplot(4,3,[1:3:7]+2);
% ax_phase_h=subplot(4,3,10+2);
plot_FRF_mag_phase(f(1:N_f_max),H_2(1:N_f_max),false,ax_mag_h,ax_phase_h,[],'\hat{H}_{2}');

set(groot,'DefaultLineLineWidth','remove');

export_figure([1:2],'==',{'CantiliverBeamImpactTest_Signals';'CantiliverBeamImpactTest_FRF'})