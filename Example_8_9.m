function Example_8_9()
clc
close all

SNR_x=1;
SNR_y=1;

D_2_1=0.2;

f_s=100*1/D_2_1;    %500 Hz
T=500*D_2_1;    %100 sec
[D_t,K_temp]=samplingParameters_T_fs(T,f_s);

rng(0);
s_vec=randn(1,K_temp);
f_c=f_s/5;  %100 Hz
[b,a]=butter(9,f_c/(f_s/2));s_vec=filtfilt(b,a,s_vec);
s_vec=s_vec-mean(s_vec); s_vec=s_vec/std(s_vec); % Makes mean(s)=0 & std(s)=1;

T_start=D_2_1;
k_start=T_start/D_t;

x_vec=s_vec(k_start+1:K_temp);
y_vec=s_vec((k_start+1:K_temp)-k_start);
K=K_temp-k_start;

rng(1);x_vec=x_vec+std(x_vec)*randn(1,K)/SNR_x;
rng(2);y_vec=y_vec+std(y_vec)*randn(1,K)/SNR_y; 

kappa_max=5*D_2_1/D_t;    %kappa_max at tau=5*D_2_1=1 sec.=T/100
r_xx=xcorr(x_vec,x_vec,kappa_max,'unbiased');
r_xy=xcorr(y_vec,x_vec,kappa_max,'unbiased');

f=(0:kappa_max-1)*f_s/kappa_max;    %Why discard the -ve tau part?
S_xy=fft(r_xy(kappa_max+1:end-1));  %Why S_xy is calculated here using fft,while in 9.3 it is calculated using cpsd?

subplot(2,1,1)
plot(f(1:kappa_max/2+1),unwrap(angle(S_xy(1:kappa_max/2+1))))
xlabel('$f$ (Hz)','interpreter','latex')
ylabel('$\angle S_{xy}(f)$ (rad)','interpreter','latex');

ind=find(f==f_c);
P1=polyfit(f(2:ind),unwrap(angle(S_xy(2:ind))),1);
hold on;plot(f(2:ind),P1(1)*f(2:ind)+P1(2));
t_delay1=-P1(1)/(2*pi)

N=2*kappa_max;
f=(0:N-1)*f_s/N;
S_xx=fft(r_xx(1:N));
% S_xx=fft(r_xx(1:N)).*exp(i*2*pi.*f*(kapa_max*D_t));   %This compensate the delay due to the inclusion of the negative part of tau

S_xy=fft(r_xy(1:N));
% S_xy=fft(r_xy(1:N)).*exp(i*2*pi.*f*(kapa_max*D_t));   %This compensate the delay due to the inclusion of the negative part of tau

H1=S_xy./S_xx;

subplot(2,1,2)
plot(f(1:kappa_max+1),unwrap(angle(H1(1:kappa_max+1))))
xlabel('$f$ (Hz)','interpreter','latex')
ylabel('$\angle H_{1}(f)$ (rad)','interpreter','latex');

P2=polyfit(f(2:ind),unwrap(angle(H1(2:ind))),1);
hold on;plot(f(2:ind),P2(1)*f(2:ind)+P2(2));
t_delay2=-P2(1)/(2*pi)
