%% ÿһƬ��������һ���������һ��ܸ˲���
clear;
load seattle_RTT;
[n1,n2,n3]=size(T);
 rate_sample=0.1:0.1:0.9; 
 [~,rate_number]=size(rate_sample);
 r=3; 
 normalize=max(T(:));
 X=T/normalize;
 Beta=0.5;
 B=zeros(n1,n2,n3);
 X_slice_rec=zeros(n1,n2,n3);
 tt=2;
%% 
for t=1:tt
    for i=1:rate_number
        m=round(0.1*i*n1*n2);
        for j=1:n3
            Xslice=X(:,:,j);
            [M,B1,Omega1]=M_Npass_sampling(m,Xslice,r,8,Beta );
            B(:,:,j)=B1;
            D2=sparse(M);
            [A2, iter2, svp2] = inexact_alm_mc(D2, 1e-4);
            X2=A2.U*A2.V';
            X2=M+X2.*(1-B1);
            X_slice_rec(:,:,j)=X2;
        end
        X_slice_rec=X_slice_rec*normalize;
        X_dif1=X_slice_rec-T;
        RSEasm(t,i)               =        norm(X_dif1(:))/(norm(T(:))); 
         % ADMM
        %M:��������    omega:����Ԫ�����     B:����Ԫ��Ϊ1,�ǳ���Ԫ��Ϊ0������
         alpha                  =        1.5;  %���ɳڲ���(alpha��ֵ����1.0��1.8֮��)
         maxItr                 =        50; % maximum iteration
         rho                    =        0.01; %�����������ղ���         
         p                      =        0.5 ;
         addpath('tSVD','proxFunctions','solvers') 
         myNorm  = 'tSVD_1' ;
         A=diag(sparse(double(B(:))));
         M1=A * X(:);        %A:n1n2n3*n1n2n3�ĳ���Ԫ��Ϊ1,�ǳ���Ԫ��Ϊ0       M1��n1n2n3*1�ĳ�������
         X1                     =        tensor_cpl_admm(A,M1,rho,alpha,[n1,n2,n3],maxItr,myNorm,0); 
         X1 = X1*normalize;
         X1                     =        reshape(X1,[n1,n2,n3])        ;
         X_dif                  =        X1-T                          ;
         RSEadmm(t,i)               =        norm(X_dif(:))/(norm(T(:)));               %relative_square_error ���ƽ�����  
     end
end
%% 
r1=zeros(1,rate_number);
r2=r1;

for k=1:tt
    r1=r1+RSEasm(k,:);
    r2=r2+RSEadmm(k,:);

end
rse1=r1/tt;
rse2=r2/tt;

 figure(1);
 plot(rate_sample,rse1,'r');
 hold on
 plot(rate_sample,rse2,'r--');
 xlabel('rate sample');
 ylabel('RSE');
 legend('alm','admm');
 
 plot(RSEasm,'r');
 hold on;
 plot(RSEadmm,'r--')