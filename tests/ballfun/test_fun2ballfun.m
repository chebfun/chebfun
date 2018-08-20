function pass = test_fun2ballfun( ) 
% Test with function (r+1)*cos(lam)*cos(th)
m = 10; n = 11; p = 12;
F = zeros(m,n,p);
F(1,floor(n/2),floor(p/2))=1/4;F(1,floor(n/2),floor(p/2)+2)=1/4;
F(1,floor(n/2)+2,floor(p/2))=1/4;F(1,floor(n/2)+2,floor(p/2)+2)=1/4;
F(2,floor(n/2),floor(p/2))=1/4;F(2,floor(n/2),floor(p/2)+2)=1/4;
F(2,floor(n/2)+2,floor(p/2))=1/4;F(2,floor(n/2)+2,floor(p/2)+2)=1/4;

f = ballfun(F);
g = ballfun(@(r,lam,th)(r+1).*cos(lam).*cos(th),[m,n,p]);

pass = isequal(f,g);
end
