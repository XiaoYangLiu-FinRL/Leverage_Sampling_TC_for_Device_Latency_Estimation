%%%% ���׶θܸ˷��������д�Լ40min
clear;
% for exp=0.05:0.05:0.3
%   T=fpower_law(6);
load seattle_RTT.mat
o=T(:,:,10:60);
clear T;
T=o;
%%%%%%%%%%%%%%%%%%T������
%   [n1,n2,n3]=size(T);
%   X1_=zeros(n1,n3,n2);
% for c=1:n2
%     X1_(:,:,n2)=T(:,n2,:);
% end
% clear T;
% T=web_click(:,:,1:50);
% V=T(1:100,1:100,:);
% clear T;
% T1=10*log10(T);
%  T1=10.^(T/10);
%  clear T;
%  T=T1;
% T=T1(1:100,1:100,:);
%%%%%%%%%%%%%%%%
r=5;
beta=0.5;               %���������ռ����
rho  =  0.01;
[n1,n2,n3]=size(T);
% meanvalue=mean(T(:));
% T=T-meanvalue;
normalize              =        max(abs(T(:)))                     ;
% normalize              =        -1   ;
D=zeros(n1,n2,n3);
F=D;
sample_rate=0.1:0.1:0.9;    %������
[~,sample_rate_number]=size(sample_rate);
RSE1=zeros(1,sample_rate_number);
RSE1sum=RSE1;
score_ct=1;                  %�ܸ˲���ѭ������
for st=1:score_ct
    for k=1:sample_rate_number
        
        X=T/normalize;
        mm=round(0.1*n1*n2*k);         %��������С����
        
        sp=randsample(n1*n2,round(mm*beta));
        B1=zeros(n1*n2,1);
        for nn=1:round(mm*beta)
            B1(sp(nn),1)=1;
        end
        B=reshape(B1,[n1,n2]);
        B2=zeros(n1,n2,n3);
        for kk=1:n3
            B2(:,:,kk)=B; %����������ֵ�����
        end
        
        X1=X.*B2;
        [UU,SS,VV]=reduce_tsvd(X1,r);
        U_h_slide={1,n1};
        V_l_slide={1,n2};
        u=zeros(1,n1);   %�洢�����ϵ��,������ֵ������Ƭ
        v=zeros(1,n2);   %�洢�����ϵ����������ֵ�Ĳ���Ƭ
        %ˮƽ��Ƭ horizon slide
        for i=1:n1
            U_h_slide{1,i}=reshape(UU(i,:,:),r,n3);
            u(1,i)=norm(U_h_slide{1,i}(:))^2;
            %      u(1,i)=norm(U_h_slide{1,i},2)^2;
        end
        %����Ƭ   lateral slide
        for j=1:n2
            V_l_slide{1,j}=reshape(VV(j,:,:),r,n3);
            v(1,j)=norm(V_l_slide{1,j}(:))^2;
            %     v(1,j)=norm(V_l_slide{1,j},2)^2;
        end
        p1=zeros(n1,n2);    %�ܸ˷���
        
        for i=1:n1
            for j=1:n2
                p1(i,j)=(u(1,i)*r/n1+v(1,j)*r/n2-(u(1,i)*r/n1)*(v(1,j)*r/n2))*(log(n1+n2)^2);
            end
        end
        p1_array=reshape(p1,[1,n1*n2]);
        [p1_descend,Idx]=sort(p1_array,'descend');
        remainer=fordifferent(Idx,sp);
        
        B3=zeros(1,n1*n2);
        for ii=1:round(mm*(1-beta))
            B3(1,remainer(ii))=1;
        end
        B4=reshape(B3,[n1,n2]);
        B5=zeros(n1,n2,n3);
        for kk=1:n3
            B5(:,:,kk)=B4; %�ܸ˲�������
        end
        BB=B2+B5;%������
        imcompleteRTT=X.*BB;
        for d=1:n3
            temp=imcompleteRTT(:,:,d);
            D(:,:,d)= NCS_vivaldi_all(temp, 8, length(temp), 32, 0);
        end
        for ii=1:n1
            for jj=1:n2
                for kk=1:n3
                    if D(ii,jj,kk)==0
                        F(ii,jj,kk)=0;
                    else
                        F(ii,jj,kk)=imcompleteRTT(ii,jj,kk)/D(ii,jj,kk);
                    end
                end
            end
        end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ADMM
        %M:��������    omega:����Ԫ�����     B:����Ԫ��Ϊ1,�ǳ���Ԫ��Ϊ0������
         alpha                  =        1.5;  %���ɳڲ���(alpha��ֵ����1.0��1.8֮��)
         maxItr                 =        50; % maximum iteration
         rho                    =        0.01; %�����������ղ���
         p                      =        0.5 ;
         addpath('tSVD','proxFunctions','solvers')
         myNorm  = 'tSVD_1' ;
         A=diag(sparse(double(BB(:))));
         M1=A * F(:);        %A:n1n2n3*n1n2n3�ĳ���Ԫ��Ϊ1,�ǳ���Ԫ��Ϊ0       M1��n1n2n3*1�ĳ�������
         X2                     =        tensor_cpl_admm(A,M1,rho,alpha,[n1,n2,n3],maxItr,myNorm,0);
         
         X2                     =        reshape(X2,[n1,n2,n3])        ;
         X3=X2.*D;
         
         X3 = X3*normalize;

         X_dif                  =        X3-T                          ;
         RSE1(1,k)               =        norm(X_dif(:))/(norm(T(:)));               %relative_square_error ���ƽ�����
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Omega=find(BB);
%         [ TC,psnr] = TCTF( X,Omega);
%         X11=TC*normalize;
%         X_dif                  =        X11-T                          ;
%         RSE1(1,k)               =        norm(X_dif(:))/(norm(T(:)));               %relative_square_error ���ƽ�����
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  if(k==5)
        %      save('tensor_denoise_200_56_recdB_0.5.mat','X1');
        %  end
    end
    RSE1sum=RSE1+RSE1sum;
end
%%%%%%%%%%%%%��ͼ
% open('Tube_Matx_score11_14.fig')
figure;
% title(['powerlaw����eΪ' num2str(exp)]);
plot(sample_rate,RSE1sum/score_ct,'r--');
%plot(sample_rate,RSE2);
xlabel('������');
ylabel('RSE:��Ծ������');



















