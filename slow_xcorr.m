function r_xy=slow_xcorr(x,y)
K=length(x);

r_xy=nan(1,2*K-1);
for kappa=-(K-1):K-1
    k_start=max(0,-kappa)+1;
    k_end=min(K-1,K-kappa-1)+1;
    r_xy(kappa+K)=sum(x(k_start:k_end).*y((k_start:k_end)+kappa));
end
