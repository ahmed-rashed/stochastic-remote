function [y_contaminated,y_pure]=contaminatedSignal(t,f0,SNR)

y_pure=sin(2*pi*f0*t);
y_contaminated=y_pure+randn(size(t))/SNR;

