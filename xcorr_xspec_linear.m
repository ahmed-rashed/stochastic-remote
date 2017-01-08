function r_xy_pad_col=xcorr_xspec_linear(x,y,scaleType) %Calculates cross-correlation and cross-Spectrum of the signals x and y

K=length(x);
K_pad=2^nextpow2(2*K-1);

x_pad_col=[x(:);zeros(K_pad-K,1)];
X_pad_col=fft(x_pad_col);
if all(x==y)    %Auto-correlation
    R_XY_pad_col=abs(X_pad_col).^2;
else            %Cross-correlation
    y_pad_col=[y(:);zeros(K_pad-K,1)];
    Y_pad_col=fft(y_pad_col);

    R_XY_pad_col=Y_pad_col.*conj(X_pad_col);
end
r_xy_pad_col=fftshift(ifft(R_XY_pad_col));

if nargin>2
    if strcmp(scaleType,'unbiased')
        k=(0:K-1).';
        k_vec=[k;ones(K_pad-2*K+1,1);k(end:-1:2)];
        r_xy_pad_col=r_xy_pad_col./(K-k_vec);
    end
else
    r_xy_pad_col=r_xy_pad_col/K;
end