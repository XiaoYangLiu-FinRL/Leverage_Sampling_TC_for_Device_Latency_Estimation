function [ Fhat ] = IS_p( F,B,t )
%IS_P �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[n1,n2]=size(F);
L=eye(n1);
xigma=1;
maxIter=2;
p=2;
for k=1:maxIter
    %%
    cvx_begin

    variable x(n1,n2) nonnegative;
%     temp=L * x;
%     S=svd(temp);
%     nuclearnorm=sum(S);
%      minimize (nuclearnorm)
    minimize ((norm( L * x,'fro')))

    subject to

        abs(x.*B-F)<=t;%�ָ�����и�Ԫ��������С��t
        
%         x>=0;

    cvx_end
    [U,S,V]=svd(x);
%     for i=1:n1
%         if S(i,i)<=xigma && S(i,i)>0
%             xigmak=S(i,i);
%             break
%         end
%             
%     end
    W=zeros(n1,n2);
    for i=1:n1
        for j=1:n2
            if i==j
            W(i,j)=(S(i,i)^p+xigma)^(-1/p);
            else
                W(i,j)=0;
            end
        end
    end
    xigma=xigma/2;
    L=U*W*U';  
end
Fhat=x;

end

