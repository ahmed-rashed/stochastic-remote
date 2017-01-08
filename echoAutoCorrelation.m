function echoAutoCorrelation()

f_tone=10;
T_tone=1/f_tone;
T_tone_total=5*T_tone;  %5 cycles tone burst

a=2;
b=1;

D_1=2*T_tone_total;
D_2=5*T_tone_total;
T_total=8*T_tone_total;

f_s=20*f_tone;
D_t=1/f_s;
N=T_total/D_t;
t=[0:N-1]*D_t;

N_tone=T_tone_total/D_t;
ToneBurst=sin(2*pi*f_tone*t(1:N_tone));

N_D_1=D_1/D_t+1;
N_D_2=D_2/D_t+1;

x=zeros(1,N);
x(N_D_1:N_D_1+N_tone-1)=a*ToneBurst;
x(N_D_2:N_D_2+N_tone-1)=x(N_D_2:N_D_2+N_tone-1)+b*ToneBurst;


rng(0);

SNR_vec=[inf,1];
N_SNR=length(SNR_vec);
figure
for n=1:N_SNR
    x_noisy=x+std(ToneBurst)/SNR_vec(n)*randn(1,N);
    subplot(2,N_SNR,2*n-1)
    plot(t,x_noisy)
    xlabel('$t$','interpreter', 'latex');
    if SNR_vec(n)==inf
        title('$x(t)=a\ s(t-\Delta_{1})+b\ s(t-\Delta_{2})$','interpreter', 'latex')
    else
        title('$x(t)=a \, s(t-\Delta_{1})+b \, s(t-\Delta_{2})+\mathrm{noise}$','interpreter', 'latex')
    end
    ylim([-4 4])

    [rxx, n_tau]=xcorr(x_noisy,x_noisy);    %autocorrelation not normalized because the signal is transient
    tau=n_tau*D_t;
    subplot(2,N_SNR,2*n)
    plot(tau,rxx)
    xlabel('$\tau$','interpreter', 'latex');
    title('$r_{xx}$','interpreter', 'latex')
    %axis([-2.5 2.5 -300 300])
end

export_figure(gcf,'==',{'EchoAutoCorrelation'})