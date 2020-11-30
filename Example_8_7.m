clc
close all
clearvars

SNR=.5;

T_echo=2;
T_burst=1;
T=T_burst+T_echo+3;
K=2^10; %K=1024
[Delta_t,f_s,D_f]=samplingParameters_T_N(T,K);
t_vec=(0:K-1)*Delta_t;
tau_vec=(-(K-1):K-1)*Delta_t;

s_fn_cvec={@(t) (chirp(t,5,1,15))};

x_vec=signalPulse(t_vec,T_burst,s_fn_cvec{1});
a=0.1;
y_vec=a*signalPulse(t_vec-T_echo,T_burst,s_fn_cvec{1});


h_vec=fliplr(x_vec); %% Matched filter

figure
subplot(2,2,1)
plot(t_vec,x_vec)
xlabel('$t$ (s)','interpreter','latex') 
ylabel('$x(t)$','interpreter','latex')

subplot(2,2,2)
plot(t_vec,h_vec)
xlabel('$t$ (s)','interpreter','latex')
ylabel('$h(t)$','interpreter','latex')

rng(0);
y_hat_vec=addNoise(y_vec,SNR);

subplot(2,2,3:4)
plot(t_vec,y_hat_vec)
xlabel('$t$ (s)','interpreter','latex') 
ylabel('$\hat{y}(t)$','interpreter','latex')

r_x_y_hat=xcorr(y_hat_vec,x_vec);

figure
subplot(2,1,1)
plot(tau_vec,r_x_y_hat)
xlabel('$\tau$ (sec.)','interpreter','latex') 
ylabel('$r_{x\hat{y}}(\tau)$','interpreter','latex')

out=conv(y_hat_vec,h_vec);
out=out(1:K); % or out=filter(h,1,y);

subplot(2,1,2)
plot(t_vec,out)
xlabel('$t$ (s)','interpreter','latex')
ylabel('Filtered signal,out(t)')