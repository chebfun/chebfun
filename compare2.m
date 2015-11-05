function compare2(f,g)
%function [t3t, t3s, t3c] = compare2(f,varargin)
% Input: an anonymous function f = @(x,y)
% Output: timings and other information for chebfun3s vs chebfun3t
pref = chebfunpref();
% if nargin>1
%     pref.eps = varargin{1};
% end

n = 20; x = linspace(-1,1,n)'; y = 0.5*x; fVals = f(x,y);
%n = 20; x = linspace(-1,1,n)'; y = x; z = x; fVals = f(x,y,z);

%fprintf('f = %s \n', func2str(f));

% tic, f2 = chebfun2(f); t2 = toc;
% if ~isempty(f2.cols)
%     f2Vals = f2(x,y);
%     err2 = norm(fVals - f2Vals)./norm(fVals);
%     rank_f2 = rank(f2);
% end
% 
% tic, f2Chopped = pivotChop(chebfun2(f)); t2Chopped = toc;
% if ~isempty(f2Chopped.cols)
%     f2ChoppedVals = f2Chopped(x,y);
%     err2Chopped = norm(fVals - f2ChoppedVals)./norm(fVals);
%     rank_f2Chopped = rank(f2Chopped);
% end
% fprintf('f2         %-4d %-2.1e\n', rank_f2, err2)
% fprintf('f2Chopped  %-4d %-2.1e\n\n', rank_f2Chopped, err2Chopped)

f2 = chebfun2(f); g2 = chebfun2(g);
tic, h2 = f2 + g2, t2 = toc
% if ~isempty(f2.cols)
%     f2Vals = f2(x,y);
%     err2 = norm(fVals - f2Vals)./norm(fVals);
%     rank_f2 = rank(f2);
% end

tic, h2New = chebfun2(@(x,y) f2(x,y) + g2(x,y)), t2New = toc
% if ~isempty(f2Chopped.cols)
%     f2ChoppedVals = f2Chopped(x,y);
%     err2Chopped = norm(fVals - f2ChoppedVals)./norm(fVals);
%     rank_f2Chopped = rank(f2Chopped);
% end
%fprintf('f2         %-4d %-2.1e\n', rank_f2, err2)
%fprintf('f2Chopped  %-4d %-2.1e\n\n', rank_f2Chopped, err2Chopped)


% compare3(@(x,y,z) sin(x + y + z));
% compare3(@(x,y,z) sin(20*x + y + z));
% compare3(@(x,y,z) sin(x + 20*y + z));
% compare3(@(x,y,z) sin(x + y + 20*z));
% compare3(@(x,y,z) sin(20*x + y + 20*z));
% compare3(@(x,y,z) sin(exp(-0.6*(x+2*y+3*z))))
% compare3(@(x,y,z) 1./(1 + x.^2 + y.^2 + z.^2))
% compare3(@(x,y,z) 1./(0.1 + x.^2 + y.^2 + z.^2))
% compare3(@(x,y,z) 1./(0.1 + x.^2 + y.^2 + 5*z.^2))
% compare3(@(x,y,z) 1./(0.01 + x.^2 + y.^2 + z.^2))
% %compare3(@(x,y,z) sin(80*x)+y+z,1e-14)
% compare3(@(x,y,z) sin(80*x)+y+z)
% %compare3(@(x,y,z) sin(x+20*y+z.^2).*exp(-(3+x+y.^2)),1e-14)
% compare3(@(x,y,z) sin(x+20*y+z.^2).*exp(-(3+x+y.^2)))
% compare3(@(x,y,z) abs(x+y+z).^11) % long !!!
% %compare3(@(x,y,z) abs(x+y+z).^11,1e-10) % long !!!
% compare3(@(x,y,z) abs(x).^11)
%compare3(@(x,y,z) besselj(0,10*sqrt(x.^2+y.^2+z.^2)),1e-12)
%compare3(@(x,y,z) besselj(0,20*sqrt(x.^2+y.^2+z.^2)),1e-10)
% compare3(@(x,y,z) sin(sqrt(x.^2 + y.^2./3 + z.^2+1)))
%compare3(@(x,y,z) tanh(x+y-.3) + cos(x.*y.*z)./(4+x-y-z),1e-13)
%compare3(@(x,y,z) x.^3.*y+sin(3*z))
%compare3(@(x,y,z) cos(6*pi*x).*exp(-pi*x.^2).*y.*z)
%compare3(@(x,y,z) exp(-100*(x.^2 + y.^2 + z.^2)))
%compare3(@(x,y,z) cosh(x+y).*sinh(z+x))
%compare3(@(x,y,z) cosh(x+y+2*z))
%compare3( @(x,y,t) 3*sin(2*x) .* sin(3*y) .* sin(t) ) % A 2D standing wave
% compare3( @(x,y,t) 3*sin(t + 2*x - 3*y ))            % A 2D travelling wave
%%
%f = @(x,y,z) sin(x + y./3 + z);
%f = @(x,y,z) cos(x+y).*sin(x + 2*y + z);
%f = @(x,y,z) sin(10*x + 5*y + 20*z);

%%compare3(@(x,y,z) besselj(0,10*sqrt(x.^2+y.^2+z.^2)))
%%compare3(@(x,y,z) besselj(0,20*sqrt(x.^2+y.^2+z.^2)))


% compare3(@(x,y,z) sin(sqrt(x.^2./25 + y.^2./4 + z.^2+1)))
%compare3(@(x,y,z) sin(exp(-0.8*(x+2*y+3*z))))
%4 seconds for chebfun3t, BUT
%110 seconds for chebfun3 (when chebfun3 did not have Phase II) but much
%faster now.

%compare3(@(x,y,z) sin(sqrt(x.^2./9 + y.^2./4 + z.^2+1)))
%compare3(@(x,y,z) sin(sqrt(x.^2./16 + y.^2./4 + z.^2+1)))

%compare3(@(x,y,z) sin(exp(-0.6*(x+2*y+3*z))))
%compare3(@(x,y,z) 1./(0.01 + x.^2 + y.^2 + z.^2))
%compare3(@(x,y,z) 1./(0.01 + x.^2 + y.^2 + z.^2),1e-14)
%compare3(@(x,y,z) 1./(0.1 + x.^2 + y.^2 + 5*z.^2))


% compare3(@(x,y,z) sin(20*x+y+z))
% compare3(@(x,y,z) sin(x+20*y+z))
% compare3(@(x,y,z) sin(x+y+20*z))

%f = @(x,y,z) cos(20*x + y + 80*z);

%compare3(@(x,y,z) sin(80*x)+y+z,1e-14)
%compare3(@(x,y,z) sin(x+20*y+z.^2).*exp(-(3+x+y.^2)),1e-14)
%compare3(@(x,y,z) abs(x+y+z).^11,1e-14) % long !!!
%compare3(@(x,y,z) abs(x).^11)

%compare3(@(x,y,z) abs(x+y+z).^11,1e-10)

%compare3(@(x,y,z) besselj(0,10*sqrt(x.^2+y.^2+z.^2)))
%compare3(@(x,y,z) besselj(0,20*sqrt(x.^2+y.^2+z.^2)))

%compare3(@(x,y,z) sin(sqrt(x.^2 + y.^2./3 + z.^2+1)))
%compare3(@(x,y,z) sin(sqrt(x.^2./9 + y.^2./4 + z.^2+1)))
%compare3(@(x,y,z) sin(sqrt(x.^2./16 + y.^2./4 + z.^2+1)))
%compare3(@(x,y,z) sin(sqrt(x.^2./25 + y.^2./4 + z.^2+1)))

%f = @(x,y,z) sin(x + y./3 + z);
%f = @(x,y,z) cos(x+y).*sin(x + 2*y + z);
%f = @(x,y,z) sin(10*x + 5*y + 20*z);