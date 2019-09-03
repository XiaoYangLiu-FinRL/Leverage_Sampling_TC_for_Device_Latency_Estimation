%% 默认取p=2
clear;
load seattle_RTT;
M=T(:,:,1);
% nmlz=max(M1(:));
% M=M1/nmlz;
[n1,n2]=size(M);
default_dimension=8;
B=zeros(n1,n2);
% Dtemp=D;
% B=D;
% for i=1:n1
%     for j=1:n2
%         D(i,j)=abs(i-j);%距离矩阵不知道对不对
%         Dtemp(i,j)=D(i,j);
%         if i==j
%             D(i,j)=0.1;
%         end
%     end
% end
beta=0.3;%采样率
% n=0;
% % for beta=0.5:0.1:0.9
% n=n+1;
m=round(n1*n2*beta);
si=randsample(n1*n2,m);
B(si)=1;
X=M.*B;
%% 计算距离矩阵%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 D = NCS_vivaldi_all(X, default_dimension, length(X), 32, 0);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


F=X./D;
F1 = fillmissing(F,'constant',0);
F2=F1;
F2(F1==inf)=0;
t=0.1;%误差容忍度
Fhat=IS_p(F2,B,t);

Mhat=Fhat.*D;
% Mhat=abs(Mhat);
M2=X+Mhat.*(1-B);
E=M2-M;
RSE=norm(E(:))/norm(M(:));
% end


% a(:,:,1)=[1 2;3 4];
% a(:,:,2)=[5 6;7 8];
% a(:,:,3)=[9 10;11 12];
% b=reshape(a,[3,2,2]);

% clear;
% fid = fopen('tm.2004-03-01.00-05-00.dat');
% data = textscan(fid,'%s');
% fclose(fid);
















