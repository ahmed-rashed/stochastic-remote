function cpsd_try()
clc
figure(5)
clf

N=1024;
x=rand(N,1);
y=rand(N,1);

X=fft(x);
Y=fft(y);

R_XY=Y.*conj(X)/N;
[Pxy,F]=cpsd(x,y,ones(size(x)),0,N, 'twosided');
n_max=length(Pxy);

f=linspace(0,1,N);
subplot(2,1,1)
semilogy(f,abs(R_XY),f(1:n_max),abs(2*pi*Pxy))
legend('mine','cpsd')

subplot(2,1,2)
plot(f,abs(R_XY(1:n_max))./abs(2*pi*Pxy))
