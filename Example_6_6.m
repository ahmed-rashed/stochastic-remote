clc
clearvars
close all

x=[1 3 5 3 1]; h=[9 7 5 3 1];

X=fft(x); H=fft(h);
yp=ifft(X.*H);
np=0:4;

figure(1)
subplot(3,1,1); stem(np, x, 'd','fill')
axis([-0.4 4.4 0 6])
xlabel('\itn'); ylabel('\itx_p\rm(\itn\rm)')

subplot(3,1,2); stem(np, h, 'fill')
axis([-0.4 4.4 0 10])
xlabel('\itn'); ylabel('\ith_p\rm(\itn\rm)')

subplot(3,1,3); stem(np, yp, 'fill')
axis([-0.4 4.4 0 90])
xlabel('\itn'); ylabel('\ity_p\rm(\itn\rm)')

Xz=fft([x zeros(1,length(h)-1)]); %% y=conv(x,h) will give the same results.
Hz=fft([h zeros(1,length(x)-1)]); %% In fact, Matlab uses this method.
yz=ifft(Xz.*Hz);                   
nz=0:8;

figure(2)
subplot(3,1,1); stem(nz, [x 0 0 0 0], 'd','fill')
axis([-0.4 8.4 0 6])
xlabel('\itn'); ylabel('\itx\rm(\itn\rm)')

subplot(3,1,2); stem(nz, [h 0 0 0 0], 'fill')
axis([-0.4 8.4 0 10])
xlabel('\itn'); ylabel('\ith\rm(\itn\rm)')

subplot(3,1,3); stem(nz, yz, 'fill')
axis([-0.4 8.4 0 90])
xlabel('\itn'); ylabel('\ity\rm(\itn\rm)')

hold on
stem(np, yp, 'fill')
legend({'Circular conv','Linear conv'})

