clc
close all
clearvars

f_s=1000;
T=5;
[Delta_t,K,Delta_f]=samplingParameters_T_fs(T,f_s);
t_vec=(0:K-1)*Delta_t;

rng(0);
s=randn(1,K);

fc=100;
[b,a] = butter(9,fc/(f_s/2));
s=filtfilt(b,a,s); 
s=s-mean(s);
s=s/std(s); % Makes mean(s)=0 & std(s)=1;

delta=0.2;
x=s(delta*f_s+1:end);
y=s(1:end-delta*f_s);

randn('state',1); nx=1*std(s)*randn(size(x));
randn('state',2); ny=1*std(s)*randn(size(y));

x=x+nx; y=y+ny; 

maxlag1=0.25*f_s; maxlag2=0.5*f_s;
[Rss,tau1]=xcorr(s,s,maxlag1,'unbiased');
[Rxy,tau2]=xcorr(y,x,maxlag2,'unbiased');
tau1=tau1*Delta_t; tau2=tau2*Delta_t;

figure(1) 
plot(tau1,Rss)
axis([-0.25 0.25 -0.4 1.2])
xlabel('Lag (\it\tau)') 
ylabel('Autocorrelation (\itR_s_s\rm(\it\tau\rm))')

figure(2)
plot(tau2(maxlag2+1:end),Rxy(maxlag2+1:end))
axis([0 0.5 -0.4 1.2])
xlabel('Lag (\it\tau)')
ylabel('Cross-correlation (\itR_x_y\rm(\it\tau\rm))')
