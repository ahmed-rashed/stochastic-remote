function R_XY_vec=ourPeriodogram(x_vec,y_vec ,win_vec)

K=length(x_vec);
if K~=numel(x_vec),error('x_vec must be a row or column vector'),end
if any(size(x_vec)~=size(y_vec)),error('x_vec and  y_vec must have the same sizes'),end

if nargin<3
    win_vec=ones(size(x_vec));
else
    if isempty(win_vec)
        win_vec=ones(size(x_vec));
    elseif any(size(x_vec)~=size(win_vec))
        error('x_vec and win_vec must have the same sizes'),
    end
end

rms_window=sqrt(mean(win_vec.^2));
x_weighted_vec=x_vec.*win_vec/rms_window;
y_weighted_vec=y_vec.*win_vec/rms_window;

X_weighted_vec=fft(x_weighted_vec);
Y_weighted_vec=fft(y_weighted_vec);
R_XY_vec=Y_weighted_vec.*conj(X_weighted_vec)/K;
