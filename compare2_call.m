clc, clear
compare2(@(x,y) cos(12+x+y))
compare2(@(x,y) cos(13+x+y))
compare2(@(x,y) cos(14+x+y))
compare2(@(x,y) sin(20+x+y))
compare2(@(x,y) abs(x+y).^11) % long !!!
compare2(@(x,y) besselj(0,10*sqrt(x.^2+y.^2)))
compare2(@(x,y) sin(sqrt(x.^2 + y.^2./3 +1)))
%%
%f = @(x,y) sin(x + y./3 + z);
%f = @(x,y) cos(x+y).*sin(x + 2*y + z);
%f = @(x,y) sin(10*x + 5*y + 20*z);

%%compare2(@(x,y) besselj(0,10*sqrt(x.^2+y.^2+z.^2)))
%%compare2(@(x,y) besselj(0,20*sqrt(x.^2+y.^2+z.^2)))


% compare2(@(x,y) sin(sqrt(x.^2./25 + y.^2./4 + z.^2+1)))
%compare2(@(x,y) sin(exp(-0.8*(x+2*y+3*z))))
%4 seconds for chebfun3t, BUT
%110 seconds for chebfun3 (when chebfun3 did not have Phase II) but much
%faster now.

%compare2(@(x,y) sin(sqrt(x.^2./9 + y.^2./4 + z.^2+1)))
%compare2(@(x,y) sin(sqrt(x.^2./16 + y.^2./4 + z.^2+1)))

%compare2(@(x,y) sin(exp(-0.6*(x+2*y+3*z))))
%compare2(@(x,y) 1./(0.01 + x.^2 + y.^2 + z.^2))
%compare2(@(x,y) 1./(0.01 + x.^2 + y.^2 + z.^2),1e-14)
%compare2(@(x,y) 1./(0.1 + x.^2 + y.^2 + 5*z.^2))


% compare2(@(x,y) sin(20*x+y+z))
% compare2(@(x,y) sin(x+20*y+z))
% compare2(@(x,y) sin(x+y+20*z))

%f = @(x,y) cos(20*x + y + 80*z);

%compare2(@(x,y) sin(80*x)+y+z,1e-14)
%compare2(@(x,y) sin(x+20*y+z.^2).*exp(-(3+x+y.^2)),1e-14)
%compare2(@(x,y) abs(x+y+z).^11,1e-14) % long !!!
%compare2(@(x,y) abs(x).^11)

%compare2(@(x,y) abs(x+y+z).^11,1e-10)

%compare2(@(x,y) besselj(0,10*sqrt(x.^2+y.^2+z.^2)))
%compare2(@(x,y) besselj(0,20*sqrt(x.^2+y.^2+z.^2)))

%compare2(@(x,y) sin(sqrt(x.^2 + y.^2./3 + z.^2+1)))
%compare2(@(x,y) sin(sqrt(x.^2./9 + y.^2./4 + z.^2+1)))
%compare2(@(x,y) sin(sqrt(x.^2./16 + y.^2./4 + z.^2+1)))
%compare2(@(x,y) sin(sqrt(x.^2./25 + y.^2./4 + z.^2+1)))

%f = @(x,y) sin(x + y./3 + z);
%f = @(x,y) cos(x+y).*sin(x + 2*y + z);
%f = @(x,y) sin(10*x + 5*y + 20*z);