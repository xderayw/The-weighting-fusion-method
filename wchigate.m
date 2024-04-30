function [gate,ntru]=wchigate(v,w,pf,nint,x0)
%计算给定虚警概率情况下，加权中心卡方分布的门限
% [gate,ntru]=wchigate(v,w,pf,nint,x0)
% 参数说明：
%  gate:        得到的门限值;
%  ntru:        所用的截断数；
%
%   v:          自由度向量
%   w:          权向量
%   pf:         虚警概率；
%   nint:       初始的截断数
%   x0:        初始门限搜索值, 0 by default

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

%如果所有的权值相同，则使用卡方分布来求恒虚警门限
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
% %去除加权中的过小量 ,是不是小于最大权的十倍
% del_indx=find(w<(max(w)/50));
% %去除占所有权和中小于90％的部分 
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
% %寻找相同权的部分
% while k<length(w)
%     i=k+1;
%    while i<=length(w)
%        if w(i)==w(k)
%            %若发现相同权者，则集中之
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
% %如果所有的权值相同，则使用卡方分布来求恒虚警门限
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
