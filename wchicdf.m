function [pf,err]=wchicdf(x,v,w,ntru,u0)
%calculate the false alarm rate of weighted chi square distribution
%   [pf,err]=wchipdf(x,v,w,ntru,u0)
%parameters:
%   x           -gate to be calculated;
%   v           -freedom vector;
%   w           -weight vector
%   ntru        -truncation number
%   u0          -see paper;
%   pf          -probability of false alarm;
%   err         -error bound

%zhou3sheng@163.com
%2007.6.27.23.00

if nargin<5
    u0=3;
end

% %若w中出现零或者很小的值，则剔除该值-cccccccccccccccccccccccccccc
zw=find(w/sum(w)<0.0001);
w(zw)=[];
v(zw)=[];

lv=length(v);

nv=sum(v);
%如果所有的系数相同，则退化为卡方分布
if all(w/w(1)==ones(size(w)))
    pf=1-chi2cdf(x/w(1),nv);
    err=0;
    return;
end

% w=w*10;
% x=x*10;

belta=(max(w)+min(w))/2;
p=nv/2+1;

if length(w)~=lv
    disp('v and w must have same length');
    gate=0;err=1;
    return;
end

u0=p/u0;
%calculate lj in the paper;
L=zeros(ntru,1);
for j=1:ntru

    for i=1:lv
        L(j)=L(j)+v(i)*((1-w(i)/belta)/(1+w(i)/belta*(p/u0-1))).^j;
    end
   L(j)=L(j)/2+(-1/(p/u0-1))^j;

end

%calculate mk in  the paper;
mk=zeros(ntru+1,1);
%calculate m0 firstly;
mk(1)=1;
for k=1:lv
    mk(1)=mk(1)*(1+w(k)/belta*(p/u0-1)).^(-v(k)/2);
end
mk(1)=mk(1)*(p/u0)^(nv/2)*2*belta*p/(p-u0);

for k=1:ntru

    for i=0:k-1
        mk(k+1)=mk(k+1)+mk(i+1)*L(k-i);
    end
    mk(k+1)=mk(k+1)/k;


end

%calculate sum of F(x) in the paper;
pf=zeros(size(x));

for k=0:ntru

    pf=pf+gamma(k+1)*mk(k+1)*gamma(nv/2+1)/gamma(nv/2+1+k)...
        *laguerre((nv+2)*x/4/belta/u0,nv/2,k);

end

pf=pf.*exp(-x/2/belta).*x.^(nv/2)/(2*belta).^(1+nv/2)/gamma(nv/2+1);

%calculate the truncation error;

if nargout>1
    eta=max(abs((1-w/belta)./(1-w/belta+w/belta*(p/u0))));

    err=0;
    for k=ntru+1:2*ntru
        err=err+eta^k*gamma(nv/2+k+1)/gamma(k+1)/gamma(nv/2+1);

    end
    err=err*exp(-x/2/belta).*x.^(nv/2)*abs(mk(1)).*...
        exp((nv+2)*x/8/belta/u0)/(2*belta)^(1+nv/2)*gamma(nv/2+1);
%    err=max(err);
%     disp('error:');
%     disp(err);
end

pf=1-pf;


return;