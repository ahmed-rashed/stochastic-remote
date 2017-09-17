function ccorr=ccorrFunc(x,y)
K=length(x);

ccorr=zeros(size(x));
for kappa=0:K-1
    for k=0:K-1
        y_ind=mod(k+kappa,K);
        
        ccorr(kappa+1)=ccorr(kappa+1)+x(k+1)*y(y_ind+1);
    end
end

