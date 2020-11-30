function Example_10_2()
clc
close all

A1=20; A2=30; f1=5; f2=15; wn1=2*pi*f1; wn2=2*pi*f2; 
zeta1=0.02; zeta2=0.01;
wd1=sqrt(1-zeta1^2)*wn1; wd2=sqrt(1-zeta2^2)*wn2;
f_s=100; T1=10; t1=[0:1/f_s:T1-1/f_s];
h=(A1/wd1)*exp(-zeta1*wn1*t1).*sin(wd1*t1) + (A2/wd2)*exp(-zeta2*wn2*t1).*sin(wd2*t1);

B1=1;
B2=0.5;
B3=0.2;
B4=0.05;
N1=fix(1.33*2/B1*f_s);
N2=fix(1.33*2/B2*f_s); 
N3=fix(1.33*2/B3*f_s);
N4=fix(1.33*2/B4*f_s);

Ns=500; Nt=N4*Ns;
rng(0); 
x=randn(1,Nt); 
y=filter(h,1,x); % we do not scale for convenience.

Gamma_1    =mscohere(x(1:Ns*N1),y(1:Ns*N1),hanning(N1),[],N4,f_s);
Gamma_2    =mscohere(x(1:Ns*N2),y(1:Ns*N2),hanning(N2),[],N4,f_s);
Gamma_3    =mscohere(x(1:Ns*N3),y(1:Ns*N3),hanning(N3),[],N4,f_s);
[Gamma_4,f]=mscohere(x(1:Ns*N4),y(1:Ns*N4),hanning(N4),[],N4,f_s);

plot(f,[Gamma_1 Gamma_2 Gamma_3 Gamma_4])
xlabel('$f$ (Hz)','interpreter','latex')
ylabel('$\gamma_{YX}^{2}(f)$','interpreter','latex')
axis([0 f_s/2 0 1.1])
