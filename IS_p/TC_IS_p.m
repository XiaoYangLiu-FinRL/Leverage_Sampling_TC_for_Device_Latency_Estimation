function [ Fhat ] = TC_IS_p( F,B,t )
%TC_IS_P 此处显示有关此函数的摘要
%   此处显示详细说明
[n1,n2,n3]=size(F);
L1=eye(n1);
L2=eye(n2);
L3=eye(n3);
xigma=1;
maxIter=2;
p=2;
for k=1:maxIter
    %%
    cvx_begin 
    
    variable x(n1,n2,n3) nonnegative ;
%     variable x(n1,n2,n3) ;
%     x1=unfold(x,1);
%     x2=unfold(x,2);
%     x3=unfold(x,3);
    minimize ( 0.3*(norm( L1 *unfold(x,1),'fro'))+0.3*(norm( L2 *unfold(x,2),'fro'))+0.4*(norm( L3 *unfold(x,3),'fro')) )
    
    subject to
    
    abs(x.*B-F)<=t;%恢复结果中各元素最大误差小于t
    

    
    cvx_end
    for l=1:3
    [U,S,~]=svd(tensor_unfolding(x,l));
    switch l
        case 1
             W=zeros(n1,n1);
             n=n1;
        case 2 
             W=zeros(n2,n2);
             n=n2;
        case 3 
             W=zeros(n3,n3);
             n=n3;
    end
   
    for i=1:n
        for j=1:n
            if i==j
                W(i,j)=(S(i,i)^p+xigma)^(-1/p);
            else
                W(i,j)=0;
            end
        end
    end
    xigma=xigma/2;
    L=U*W*U';
    switch l
        case 1
            L1=L;
        case 2 
            L2=L;
        case 3 
            L3=L;
    end
    end
Fhat=x;

end

