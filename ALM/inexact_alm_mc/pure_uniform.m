clc;
clear;
%������������
r=5; % %�����ȣ��ѵ�rҪ��֪��������ԭʼ��û��ȡ��ʱ���û����������衣
n1=200;
n2=200;            
error_bound=0.01; %һ�������У������С�ڶ��٣���Ϊ�ָ���ȷ
success_ratio=0.9; %ĳһԪ�أ������������У��ָ���ȷ�Ĵ���ռ�ܴ����ı�������ĳһֵʱ����Ϊ���Ԫ�ػָ���ȷ
Beta=0.9; 
alpha=0.7;
tt=2;  %���ж�����ֵ

% %�����ʸı�
 rate_sample=0.1:0.15:0.9;                  %������
 [q,rate_number]=size(rate_sample);
% rate_number=9;

%�����ʸı�
% rate_sample=[25 30 32 35 40];                      %������
% [q,rate_number]=size(rate_sample); 
%  a=n1*log10(n1);
 
  %����ԭʼ���� 

U_original=normrnd(0,1,n1,r); %������̬�ֲ������
V_original=normrnd(0,1,n2,r);
D=zeros(n1,n2);
ratio_r=zeros(n1,n2);
for i=1:n1
    D(i,i)=1/(i^alpha);
end
 X=D*U_original*V_original'*D;

% X=U_original*V_original';
sum_X=norm(X,'fro');


%Ԥ����洢���������
sum_E2=zeros(1,rate_number);
relative_error_ratio_uniform=zeros(1,rate_number);

E_all=zeros(n1,n2);
E_mean=zeros(n1,n2);
success_number_ratio=zeros(1,rate_number);



success=cell(1,tt);
exact_recovery =cell(1,tt);
all_exact=cell(1,tt);

final_success=zeros(n1,n2);

for i=1:tt
    success{1,i}=zeros(n1,n2);
    exact_recovery{1,i}=zeros(n1,n2);
    all_exact{1,i}=ones(n1,n2);
end

   for l=1:rate_number
     m=floor(rate_sample(l)*n1*n2); 
%     m=floor(rate_sample(l)*a);       %�ܵ�����     %�ܵ�������
    final_success=zeros(n1,n2);          %ÿ�ζ�Ҫ��ʼ��һ�¡�
    
    
    E_all=zeros(n1,n2);                  % E_all����ÿһ�������ʶ�Ӧ��tt������Ĵ����ܺ�
    
    
    for i=1:tt
    
    %�۲⼯��
     [M2,Omega2]=M_uniform(m,X);
%    [M2,omega2] = my_randsample_unif(m,X);
    
    %ָʾ�� 
 %   B2=B_product(n1,n2,Omega2);
   
   
    D2=sparse(M2);
    

    [A2 iter2 svp2] = inexact_alm_mc(D2, 1e-4);
   
 
    X2=A2.U*A2.V';
    
    %�����
    E2=X-X2;
    
    %������ȷ������������С��0.05Ϊ��ȷ��׼��
    success{1,i}=abs(E2)./abs(X);
    true_position=find(success{1,i}<=error_bound); %�������С��5%������ȷ
    exact_recovery{1,i}(true_position)=all_exact{1,i}(true_position);
    
    %������ȷ���ܴ���
    final_success= final_success+exact_recovery{1,i};
    
    %����tt�ε������
    E_all=E_all+E2;
    
    
    end
    
    % ����ȷ����Ϊ������׼
     really_success=find(final_success>=success_ratio*tt); %10������7�γɹ�����ɹ�
     success_number=size(really_success);  %���������ɹ���Ԫ�ظ���
     success_number_ratio(1,l)=success_number(1,1)/(n1*n2);
     
    % ��������Դ���Ϊ������׼
     E_mean=E_all/tt; % E_mean����ÿһ�������ʶ�Ӧ��tt������Ĵ����ֵ
     sum_E2(1,l)=norm(E_mean,'fro');
     relative_error_ratio_uniform(1,l)=sum_E2(1,l)/sum_X; %���������ռԭʼ��ķ�����

    
  end


%��ͼ
%�����ʺ���ȷ�ָ������Ĺ�ϵ
 figure(1);
 plot(rate_sample,success_number_ratio,'r-*');
 xlabel('rate_sample');
 ylabel('��ȷ�ָ�����');
 
 %�����ʺ���Դ����ʵĹ�ϵ
 figure(2)
 plot(rate_sample,relative_error_ratio_uniform,'b-*');
 xlabel('rate_sample');
 ylabel('��Դ������');
