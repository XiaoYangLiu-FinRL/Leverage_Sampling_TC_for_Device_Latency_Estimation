clc;
clear;
%基本参数设置
r=5; % %矩阵秩，难道r要已知。。对于原始阵还没获取的时候就没法做这个假设。
n1=200;
n2=200;            
error_bound=0.01; %一次试验中，当误差小于多少，认为恢复正确
success_ratio=0.9; %某一元素，在所有试验中，恢复正确的次数占总次数的比例大于某一值时，认为这个元素恢复正确
Beta=0.9; 
alpha=0.7;
tt=2;  %运行多次求均值

% %采样率改变
 rate_sample=0.1:0.15:0.9;                  %采样率
 [q,rate_number]=size(rate_sample);
% rate_number=9;

%采样率改变
% rate_sample=[25 30 32 35 40];                      %采样率
% [q,rate_number]=size(rate_sample); 
%  a=n1*log10(n1);
 
  %产生原始数据 

U_original=normrnd(0,1,n1,r); %生成正态分布随机数
V_original=normrnd(0,1,n2,r);
D=zeros(n1,n2);
ratio_r=zeros(n1,n2);
for i=1:n1
    D(i,i)=1/(i^alpha);
end
 X=D*U_original*V_original'*D;

% X=U_original*V_original';
sum_X=norm(X,'fro');


%预分配存储结果的向量
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
%     m=floor(rate_sample(l)*a);       %总的样本     %总的样本数
    final_success=zeros(n1,n2);          %每次都要初始化一下。
    
    
    E_all=zeros(n1,n2);                  % E_all代表每一个采样率对应的tt次试验的错误总和
    
    
    for i=1:tt
    
    %观测集合
     [M2,Omega2]=M_uniform(m,X);
%    [M2,omega2] = my_randsample_unif(m,X);
    
    %指示阵 
 %   B2=B_product(n1,n2,Omega2);
   
   
    D2=sparse(M2);
    

    [A2 iter2 svp2] = inexact_alm_mc(D2, 1e-4);
   
 
    X2=A2.U*A2.V';
    
    %误差项
    E2=X-X2;
    
    %计算正确个数比例（以小于0.05为正确标准）
    success{1,i}=abs(E2)./abs(X);
    true_position=find(success{1,i}<=error_bound); %错误误差小于5%就算正确
    exact_recovery{1,i}(true_position)=all_exact{1,i}(true_position);
    
    %计算正确的总次数
    final_success= final_success+exact_recovery{1,i};
    
    %计算tt次的总误差
    E_all=E_all+E2;
    
    
    end
    
    % 以正确个数为衡量标准
     really_success=find(final_success>=success_ratio*tt); %10次里有7次成功就算成功
     success_number=size(really_success);  %计算真正成功的元素个数
     success_number_ratio(1,l)=success_number(1,1)/(n1*n2);
     
    % 以总体相对错误为衡量标准
     E_mean=E_all/tt; % E_mean代表每一个采样率对应的tt次试验的错误均值
     sum_E2(1,l)=norm(E_mean,'fro');
     relative_error_ratio_uniform(1,l)=sum_E2(1,l)/sum_X; %计算误差阵占原始阵的范数比

    
  end


%画图
%采样率和正确恢复比例的关系
 figure(1);
 plot(rate_sample,success_number_ratio,'r-*');
 xlabel('rate_sample');
 ylabel('正确恢复比例');
 
 %采样率和相对错误率的关系
 figure(2)
 plot(rate_sample,relative_error_ratio_uniform,'b-*');
 xlabel('rate_sample');
 ylabel('相对错误比例');
