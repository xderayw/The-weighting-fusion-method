function y=laguerre(x,alpha,j)
%calculate the generalized laguerre poynomial 

%2007.7.28.09.16

y=zeros(size(x));

for k=0:j
    
    y=y+gamma(j+alpha+1)*(-x).^k/gamma(j-k+1)/gamma(alpha+k+1)/gamma(k+1);
end

