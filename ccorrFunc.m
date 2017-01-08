function ccorrRes = ccorrFunc(x,y)
K=length(x);

ccorrRes=zeros(size(x));
for kappa=0:K-1
    for k=0:K-1
        y_ind=mod(k+kappa,K);
        
        ccorrRes(kappa+1)=ccorrRes(kappa+1)+x(k+1)*y(y_ind+1);
    end
end

