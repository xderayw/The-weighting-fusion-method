function [gate,ntru]=wchigate(v,w,pf,nint,x0)
%��������龯��������£���Ȩ���Ŀ����ֲ�������
% [gate,ntru]=wchigate(v,w,pf,nint,x0)
% ����˵����
%  gate:        �õ�������ֵ;
%  ntru:        ���õĽض�����
%
%   v:          ���ɶ�����
%   w:          Ȩ����
%   pf:         �龯���ʣ�
%   nint:       ��ʼ�Ľض���
%   x0:        ��ʼ��������ֵ, 0 by default

if nargin<5
    x0=0;
end

if pf>=1
    disp('pf must a value less than 1');
    gate=0;ntru=1;
    return;
end
%test the length

nv=length(v);
if length(w)~=nv
    disp('v and w must have same length');
    gate=-1 
    err=1;
    return;
end

%������е�Ȩֵ��ͬ����ʹ�ÿ����ֲ�������龯����
if all(w/max(w)==ones(size(w)))
    gate=chi2inv(1-pf,sum(v))*w(1);
%    err=0;
    return;
end
step=1;
for k=1:4
    step=step/10;
    x0=wchigatepf(x0,v,w,nint,pf,step);
end
gate=x0;
return;

xpf=wchicdf(x0,v,w,nint)-pf;
if xpf>0
    while xpf>0
        x0=x0+step;
        xpf=wchicdf(x0,v,w,nint)-pf;
    end
elseif xpf<0
    while xpf<0
        x0=x0-step;
        xpf=wchicdf(x0,v,w,nint)-pf;
    end
else
    gate=x0;
    err=0;
    return;
end
return;
% 
% %ȥ����Ȩ�еĹ�С�� ,�ǲ���С�����Ȩ��ʮ��
% del_indx=find(w<(max(w)/50));
% %ȥ��ռ����Ȩ����С��90���Ĳ��� 
% % k=1;sumW=sum(w);sortW=sort(w);
% % while(sum(sortW(1:k))<sumW*.9)
% %     k=k+1;
% % end
% % del_indx=find(w<=sortW(k));
% %delete;
% w(del_indx)=[];
% v(del_indx)=[];
% 
% k=1;
% %Ѱ����ͬȨ�Ĳ���
% while k<length(w)
%     i=k+1;
%    while i<=length(w)
%        if w(i)==w(k)
%            %��������ͬȨ�ߣ�����֮
%            w(i)=[];
%            v(k)=v(k)+v(i);
%            v(i)=[];
%            nv=nv-1;
%            i=i-1;
%        end
%        i=i+1;
%    end
%     k=k+1;
% end
% nv=length(v);
% %������е�Ȩֵ��ͬ����ʹ�ÿ����ֲ�������龯����
% if all(w/w(1)==ones(size(w)))
%     gate=chi2inv(1-pf,v)*w(1);
%     err=0;
%     return;
% end
% %initialize the search procedure;
% if nargin<5
%     x1=0;
% %     while wchicdf(x1,v,w,nint)>pf
% %         x1=x1+0.1;
% %     end
% %     x0=x1-0.1;
% else
%     x1=x0;
%     
% %     if wchicdf(x1,v,w,nint)>pf
% %         while wchicdf(x1,v,w,nint)>pf
% %             x1=x1+.1;
% %         end
% %         x0=x0-.1;
% %     elseif wchicdf(x1,v,w,nint)<pf
% %         while wchicdf(x1,v,w,nint)<pf
% %             x1=x1-.1;
% %         end
% %         x0=x1;x1=x0+.1;
% %     else
% %         g=x1;
% %         return;
% %     end
% end
% % ntru=nint;
% % while wchierr(x1,v,w,ntru)>pf/2
% %     ntru=ntru+1;
% % end
% % gate=fzero(@(x) wchicdf(x,v,w,ntru)-pf,x1,optimset('TolFun',pf/20));
% % return;
% % 
% odpf=5;
% xpf=wchicdf(x1,v,w,nint)-pf;
% if xpf>0
%     while xpf>0
%         x1=x1+1;
%         xpf=wchicdf(x1,v,w,nint)-pf;
% %         if xpf<0
% %             while wchierr(x1,v,w,ntru)>pf/2
% %                 ntru=ntru+1;
% %             end
% %         end
%     end
%     x0=x1-1;
%     
% elseif xpf<0
%     while xpf<0
%         x1=x1-1;
% 
%         xpf=wchicdf(x1,v,w,nint)-pf;
% %         if xpf>0
% %             while(wchierr(x1,v,w,ntru)>pf/2)
% %                 ntru=ntur+1;
% %             end
% %         end
% 
%     end
%     x0=x1;
%     x1=x1+1;
% else
%     if wchierr(x1,v,w,nint)<pf/odpf
%         gate=x1;
%         ntru=nint;
%         return;
%     end
% end
% 
% xpf=1;
% while abs(xpf)>pf/odpf
%     %calculate the order or truncation number of the target ;
%     x2=(x0+x1)/2;
% 
%     xpf=wchicdf(x2,v,w,nint)-pf;
% 
%     if xpf>0
%         x0=x2;
%     elseif xpf<0
%         x1=x2;
%     else
%         break;
%     end
% end
% 
% gate=x2;
% 
% ntru=nint;
% return;
%EOF#
