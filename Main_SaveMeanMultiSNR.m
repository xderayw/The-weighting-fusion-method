close all;clear;
%保存多个展开点
%均值eta01
tic
pf=1e-4;
K=[32 29 26 23];%reference number
%K=[32 32 32 32];%reference number
L=[16 16 16 16];%snapshot
site_num=length(K);
M=K+1-L;
rou=(K+2-L)./(K+1);%损失因子均值
w_modify=(K+1)./(K+1-L);
trial_num0=1e6;
trial_num1=2e6;
SNR=0:1:20;
lamda_snr=10.^(SNR/10);
snr_num=length(lamda_snr);

lamda_snrtimes=[1 2 1 5]; %局部信噪比的比例
lamda_snrtimes=[10 6 2 1]; %局部信噪比的比例
lamda_snrtimes=[1 7 1 1]; %局部信噪比的比例
 lamda_snrtimes=[1 1 14 15]; %局部信噪比的比例

snr_matrix=[];
for i=1:length(lamda_snrtimes)
    lamda_snr_local=lamda_snrtimes(i)*lamda_snr/sum(lamda_snrtimes);
    beta_snr_local=lamda_snr_local/rou(i);
    snr_matrix=[snr_matrix beta_snr_local'];
end
snr_choice=1:snr_num;
%snr_choice=9;
snr_choice_num=length(snr_choice);

r1=zeros(1,trial_num1);
r_expand1=zeros(site_num,snr_choice_num);
eta01=zeros(site_num,snr_choice_num);
for j=1:snr_choice_num
    for i=1:site_num
        randsig=(randn(M(i)+1,trial_num1)+1i*randn(M(i)+1,trial_num1))/sqrt(2);
        loss_factor=betarnd(K(i)-L(i)+2,L(i)-1,1,trial_num1);
        r1=(K(i)+1)*log(1+abs(randsig(1,:)+sqrt(loss_factor)*sqrt(snr_matrix(j,i))).^2./sum(abs(randsig(2:end,:)).^2));
        %r0=(K(i)+1)*log(1+abs(randsig(1,:)).^2./sum(abs(randsig(2:end,:)).^2));
        r_expand1(i,j)=mean(r1);
        eta01(i,j)=1-exp(-r_expand1(i,j)./(K(i)+1));
    end
end
 % file_name=['MultiMeanValue1215' '.mat'];
% file_name=['MultiMeanValue10621' '.mat'];
% file_name=['MultiMeanValue1711' '.mat'];
 file_name=['MultiMeanValue111415' '.mat'];

save(file_name,'eta01','snr_matrix')
toc