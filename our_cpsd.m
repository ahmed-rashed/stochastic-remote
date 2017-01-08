function R_XY_vec=our_cpsd(x_vec, y_vec ,window ,K_overlap)

K_total=length(x_vec);
K=length(window);
x_size=size(x_vec);
ind_1=find(x_size==1,1);

if isempty(ind_1),error('x_vec must be a row or column vector');end
if size(y_vec,ind_1)~=1,error('x_vec and  y_vec must be row or column vectors'),end
if length(y_vec)~=K_total,error('x_vec and y_vec must have the same size'),end
if size(window,ind_1)~=1,error('window and x_vec must be row or column vectors'),end
if K>K_total,error('window length must be smaller than or equals to x_vec length'),end
if K_overlap>=K,error('N_overlap must be smaller than the window length'),end

Q=floor((K_total-K_overlap)/(K-K_overlap));
R_XY_vec=zeros(size(window));
for q=1:Q
    k_start=(q-1)*(K-K_overlap)+1;
    R_XY_vec=R_XY_vec+ourPeriodogram(x_vec(k_start:k_start+K-1),y_vec(k_start:k_start+K-1),window);
end
R_XY_vec=R_XY_vec/Q;