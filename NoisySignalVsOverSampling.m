function NoiseVsOverSampling()
clc

T=1;
df=1/T;

f0=20*1/T;

SNR=.5;
N_vec=2.^[7:9];
figure(1);clf
for n=1:length(N_vec)
    N=N_vec(n);

    dt=T/N;
    f_s=N*df;

    n_f_max=floor(N/2);
    f_max=n_f_max*df;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t=[0:N-1]*dt;
    f=[0:N-1]*df;

    [x,x_pure]=contaminatedSignal(t,f0,SNR);
    X=fft(x);

    subplot(length(N_vec),2,(n-1)*2+1)
    plot(t,x)
    hold on;
    plot(t,x_pure,'k','LineWidth',2);
    xlabel('$t$','interpreter','latex')
    title(['$x(t)=\sin(2\pi f_0 t)+n_{\mathrm{white}},:f_0=',num2str(f0),'$ Hz,SNR=',num2str(SNR),' \& $N=',int2str(N),'$'],'interpreter','latex')
    legend({'sin + noise','pure sin'});
    ylim([min(x),max(x)])

    subplot(length(N_vec),2,(n-1)*2+2)
    plot(f(1:n_f_max+1),abs(X(1:n_f_max+1))*2/N);
    title(['$|X(f)| \times 2 / N$,SNR=',num2str(SNR),' \& $N=',int2str(N),'$'],'interpreter','latex')
    xlabel('$f$','interpreter','latex')
    xlim([0,f_max])
    ylim([0,1.5])
    set(gca,'YGrid','on')
    hold on;
    plot([f0,f0],get(gca,'YLim'),'--m')
    hold off
end
