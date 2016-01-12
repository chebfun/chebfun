function pass = test_roots10( pref ) 
% Check that the marching squares and Bezoutian agree with each other. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.eps; 
j = 1;

%% (Marching Squares misses some solution along edge)
f = chebfun2(@(x,y)(x-1).*(cos(x.*y.^2)+2)); 
g = chebfun2(@(x,y)sin(8*pi*y).*(cos(x.*y)+2)); 
r2 = roots([f;g],'resultant'); 
pass(j) = ( norm(sort(r2(:,1))-1) < tol ); j = j + 1; 
pass(j) = ( norm(sort(r2(:,2))-linspace(-1,1,17).') < tol ); j = j + 1;

%%  So many good test problems.  Put into test later.
% p1=@(x,y) (x.^2-1).^2-y.^2.*(3+2.*y); %pretzel
% rr = randn(1,2); p2=@(x,y) y.^2-x.^3+5*rr(1).*x+5*rr(2); %elliptic
% p3=@(x,y) prod(randn(d,1).*x-randn(d,1).*y); %star; origin is highly singular and may cause problems
% p4=@(x,y) (x.^2+y.^2+x.*y).*[y.^2 y 1].*randn(3).*[x.^2;x;1]+[y.^2 y 1].*100.*eps*randn(3).*[x.^2;x;1]; %perturbation of a given curve; two copies will be close to having a common component
% p5=@(x,y) 1e-7*(((3.*x+3.*y).^5-54.*x.*y).*((3.*x-0.25).^6+(3.*y+0.36).^7-0.1.*(3.*x+3.*y).^8).*(x-y).*(x+y).*x.*y.*((3.*x+3.*y).^4+(3.*x+3.*y).^3+(3.*x-3.*y).^9).*(x-0.5).*(x-1).*(x+0.5).*(x+1).*(y-0.5).*(y-1).*(y+0.5).*(y+1).*(x-y+0.5).*(x-y-0.5).*(x+y-0.5).*(x+y+0.5)); %embarassingly full of singular points
% p5=@(x,y) 1e-4*(((3.*x+3.*y).^5-54.*x.*y).*(x-y).*(x+y).*x.*y.*((3.*x+3.*y).^4+(3.*x+3.*y).^3+(3.*x-3.*y).^9).*(x-0.5).*(x-1).*(x+0.5).*(x+1).*(y-0.5).*(y-1).*(y+0.5).*(y+1).*(x-y+0.5).*(x-y-0.5).*(x+y-0.5).*(x+y+0.5)); %embarassingly full of singular points
% cp=-cos([0:d]*pi/d);p6=@(x,y) prod(x-cp).*prod(y-cp); %grid such that each chebyshev point is singular
% a=1/4;p7=@(x,y) (x.^2 + y.^2 - 2.*a.*x).^2  - 4*a^2.*(x.^2 + y.^2); %cardioid
% p8=@(x,y) (y.^2-x.^2).*(x-1).*(2.*x-3)-4.*(x.^2+y.^2-2.*x).^2; %Ampersand
% p9=@(x,y) (x.^2-1/9).*(x-1/3).^2+(y.^2-1/9).^2; %two cusps
% p10=@(x,y) x.^4+2.*x.^2.*y.^2+y.^4-x.^3-3.*x.*y.^2; %Irish curve
% p11=@(x,y) (x.^2+y.^2).^3-4.*x.^2.*y.^2; %lucky Irish curve
% p12=@(x,y) (-6.*x+9.*x.^2+9.*y.^2).^2-9.*x.^2-9.*y.^2; %trisectrix
% p13=@(x,y) x.^6+y.^6-x.^2; %butterfly
% p14=@(x,y) x.^3.*y+y.^3+x; %Klein curve
% r = rand; 
% p15=@(x,y) ((x-0.5).^2+y.^2).*((x+0.5).^2+y.^2)-r; %Cassini oval
% p16=@(x,y) (16*x.^2+16*y.^2-4).^3-1728.*y.^2; %nephroid
% p17=@(x,y) (x.^2+y.^2).^2-2.*x.^2-2.*y.^2; %Bernoulli lemniscate
% p18=@(x,y) (16*(x.^2)+16*(y.^2)).^2+288.*(x.^2+y.^2)-8.*(64*(x.^3)-192*x.*(y.^2))-27; %deltoid
% p19=@(x,y) (y.^2).*(y.^2-0.5)-x.^2.*(x.^2-1); %devil's curve
% a=0:1/10:1; 
% p20a=@(x,y) prod(4.*y.^2.*(4.*y.^2-a.^2)-4.*x.^2.*(4.*x.^2-1)); %diabolic combination of devil's curves
% p20b=@(x,y) (4*(y.^2).*(4*(y.^2)-a(1))-4.*(x.^2).*(4*(x.^2)-1))... 
% .*(4*(y.^2).*(4*(y.^2)-a(2))-4.*(x.^2).*(4*(x.^2)-1))...
% .*(4*(y.^2).*(4*(y.^2)-a(3))-4.*(x.^2).*(4*(x.^2)-1))...
% .*(4*(y.^2).*(4*(y.^2)-a(end-1))-4.*(x.^2).*(4*(x.^2)-1))...
% .*(4*(y.^2).*(4*(y.^2)-a(end))-4.*(x.^2).*(4*(x.^2)-1));
% p20=@(x,y) (4*(y.^2).*(4*(y.^2)-a(1))-4.*(x.^2).*(4*(x.^2)-1))... 
% .*(4*(y.^2).*(4*(y.^2)-a(end))-4.*(x.^2).*(4*(x.^2)-1));
% 
% f = @(x,y) p1(x,y);g = @(x,y) p2(x,y); 
% f = @(x,y) p5(x,y);g = @(x,y) p7(x,y); numsol = 30;
% f = @(x,y) p20(x,y);g = @(x,y) p18(x,y); numsol = 12; %keep
% f = @(x,y) p8(x,y);g = @(x,y) p9(x,y); %double root
% f = @(x,y) p10(x,y);g = @(x,y) p12(x,y);  % double at boundary x= 1
% f = @(x,y) p11(x,y);g = @(x,y) p13(x,y);  % highly multiple at 0
% f = @(x,y) p14(x,y);g = @(x,y) p15(x,y);  % fine
% f = @(x,y) p16(x,y);g = @(x,y) p17(x,y);  % none
% f = @(x,y) p18(x,y);g = @(x,y) p19(x,y);  % cool

end
