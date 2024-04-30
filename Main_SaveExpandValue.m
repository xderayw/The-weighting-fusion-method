close all;clear;
%������չ����
%��ֵ�Լ�H1�·�λ��
%H0�·�λ��
tic
pf=1e-4;
%
K=[32 29 26 23];%reference number
L=[16 16 16 16];%snapshot
lamda_snrtimes=[1 7 6 2]; %�ֲ�����ȵı���
% %
% K=[32 29 26 23];%reference number
% L=[16 16 16 16];%snapshot
% lamda_snrtimes=[10 6 2 1]; %�ֲ�����ȵı���
% %
% K=[32 29 26 23];%reference number
% L=[16 16 16 16];%snapshot
% lamda_snrtimes=[1 1 14 15]; %�ֲ�����ȵı���
% %
% K=[32 32 32 32];%reference number
% L=[16 16 16 16];%snapshot
% lamda_snrtimes=[1 3 5 7]; %�ֲ�����ȵı���
% %
% K=[32 32 32 32];%reference number
% L=[8 10 16 20];%snapshot
% lamda_snrtimes=[8 4 2 1]; %�ֲ�����ȵı���
% %
% K=[32 23];%reference number
% L=[16 16];%snapshot
% lamda_snrtimes=[8 1]; %�ֲ�����ȵı���



site_num=length(K);
M=K+1-L;
rou=(K+2-L)./(K+1);%��ʧ���Ӿ�ֵ
w_modify=(K+1)./(K+1-L);
trial_num0=1e6;
trial_num1=5e6;
SNR=0:1:20;
lamda_snr=10.^(SNR/10);
snr_num=length(lamda_snr);


snr_matrix=[];
for i=1:length(lamda_snrtimes)
    lamda_snr_local=lamda_snrtimes(i)*lamda_snr/sum(lamda_snrtimes);
    beta_snr_local=lamda_snr_local/rou(i);
    snr_matrix=[snr_matrix beta_snr_local'];
end
snr_matrix=snr_matrix(10,:);%dB

r_expand1=zeros(1,site_num);
r_expand2=zeros(1,site_num);
r_expand3=zeros(1,site_num);
r1=zeros(1,trial_num1);
pro_standard=[0.9];
%pro_standard=[0.9 0.7 0.5 0.3];%H1�·�Ϊ��
expand_num=length(pro_standard);

eta02=zeros(expand_num,site_num);
eta03=zeros(expand_num,site_num);
eta04=zeros(expand_num,site_num);

for j=1:expand_num
    for i=1:site_num
        randsig=(randn(M(i)+1,trial_num1)+1i*randn(M(i)+1,trial_num1))/sqrt(2);
        loss_factor=betarnd(K(i)-L(i)+2,L(i)-1,1,trial_num1);
        r1=(K(i)+1)*log(1+abs(randsig(1,:)+sqrt(loss_factor)*sqrt(snr_matrix(i))).^2./sum(abs(randsig(2:end,:)).^2));
        r0=(K(i)+1)*log(1+abs(randsig(1,:)).^2./sum(abs(randsig(2:end,:)).^2));
        r_expand1(i)=mean(r1);
        rs1=sort(r1);
        rs0=sort(r0);
        r_expand2(i)=rs1(trial_num1-trial_num1*pro_standard(j));
        r_expand3(i)=rs0(trial_num1-trial_num1*pro_standard(j));
        r_expand4(i)=(rs1(trial_num1/2)+rs1(trial_num1/2+1))/2;
    end
    %%չ�����ѡȡ 
    eta02(j,:)=1-exp(-r_expand2./(K+1));%��λ��չ����H1
    eta03(j,:)=1-exp(-r_expand3./(K+1));%��λ��չ����H0
    eta04(j,:)=1-exp(-r_expand4./(K+1));
end
eta01=1-exp(-r_expand1./(K+1));%��ֵ--չ����
file_name=['ExpandValue1762' '.mat'];
% file_name=['ExpandValue10621' '.mat'];
% file_name=['ExpandValue1145' '.mat'];
% file_name=['ExpandValue1357' '.mat'];
%file_name=['ExpandValue8421' '.mat'];
% file_name=['ExpandValue81' '.mat'];


save(file_name,'eta01','eta02','eta03','eta04')
toc