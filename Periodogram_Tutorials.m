clc
close all
clearvars
set(groot,'DefaultLineLineWidth',1);

rng(0);
f_s = 1000;
N=500;
K=N;
D_t=1/f_s;
T=N*D_t;
D_f=f_s/N;
t_vec = [0:K-1]*D_t;
f_vec = [0:N-1]*D_f;
A_signals=[1 2];
SNR=10;

f_0_signals=(floor([150;350]/D_f))*D_f; %This is to avoid leakage
f_0_signals(2)=f_0_signals(2)+D_f/2;    %This is to obtain leakage in the 2nd frequency

x_vec=A_signals*sin(2*pi*f_0_signals*t_vec);
x_vec=x_vec+std(x_vec)/SNR*randn(1,N);

windows=[rectwin(N).';hann(N).'];
windowType={'','Hann'};
for ii=1:size(windows,1)
    figure
    
    %Our Periodogram implementation
    R_XX_our=ourPeriodogram(x_vec, x_vec ,windows(ii,:));
    subplot(2,1,1)
    semilogy(f_vec(1:floor(N/2)),R_XX_our(1:floor(N/2))); grid on;
    xlabel('$f$ (Hz)', 'interpreter', 'latex');
    if isempty(windowType{ii})
        ylabel('$R_{XX}$ (dB)', 'interpreter', 'latex');
        titleString='Our periodogram implementation, $R_{XX}$';
    else
        ylabel('$R_{\widetilde{X}\widetilde{X}}$ (dB)', 'interpreter', 'latex');
        titleString=['Our modified periodogram implementation, $R_{\widetilde{X}\widetilde{X}}$, using ',windowType{ii},' window'];
    end
    title(titleString, 'interpreter', 'latex');

    %Matlab Periodogram implementation
    R_XX_Matlab = periodogram(x_vec,windows(ii,:),N,f_s,'twosided')*f_s;
    subplot(2,1,2)
    semilogy(f_vec(1:floor(N/2)),R_XX_Matlab(1:floor(N/2))); grid on;
    xlabel('$f$ (Hz)', 'interpreter', 'latex');
    if isempty(windowType{ii})
        ylabel('$R_{XX}$ (dB)', 'interpreter', 'latex');
        titleString='Matlab periodogram implementation, $R_{XX}$';
    else
        ylabel('$R_{\widetilde{X}\widetilde{X}}$ (dB)', 'interpreter', 'latex');
        titleString=['Matlab modified periodogram implementation, $R_{\widetilde{X}\widetilde{X}}$ using ',windowType{ii},' window',];
    end
    title(titleString, 'interpreter', 'latex');

    
    Error=max(abs(R_XX_our'-R_XX_Matlab))
end

export_figure([1:2],'',{'PeriodogramTutorial';'ModifiedPeriodogramTutorial'});

set(groot,'DefaultLineLineWidth','remove');