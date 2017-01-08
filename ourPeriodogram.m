function R_XY_vec=ourPeriodogram(x_vec, y_vec ,window)

K=length(x_vec);

if K~=numel(x_vec),error('x_vec must be a row or column vector'),end
if any(size(x_vec)~=size(y_vec)),error('x_vec and  y_vec must have the same sizes'),end
if any(size(x_vec)~=size(window)),error('x_vec and  window must have the same sizes'),end

rms_window=sqrt(mean(window.^2));
x_weighted_vec=x_vec.*window/rms_window;
y_weighted_vec=y_vec.*window/rms_window;

X_weighted_vec = fft(x_weighted_vec);
Y_weighted_vec = fft(y_weighted_vec);
R_XY_vec = Y_weighted_vec.*conj(X_weighted_vec)/K;
