function pass = test_roots_syntax( pref )
% Check chebfun2v rootfinding syntax. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.eps; 
j = 1;

% This is James Whidborne's problem that caused a silly error with the
% rootfinding syntax in Chebfun: 
% set up to solve a common roots problem
pdomain = [0.01 25]; % p domain eq31
phidomain = [0 2*pi]; % phi domain

V1= [phidomain  pdomain];

phi=chebfun2(@(phi,p) phi, V1);
p=chebfun2(@(phi,p) p, V1);

q=1.4; %

A=(q*cos(q*pi).*cos(phi)+(2*q^2-1)*sin(q*pi).*sin(phi))/(q^2-1)...
        + pi*sin(q*pi)./p + cos(q*pi)./(q*p); %eq32
B=(-q*sin(q*pi).*cos(phi)+(2*q^2-1)*cos(q*pi).*sin(phi))/(q^2-1)...
        + pi*cos(q*pi)./p -sin(q*pi)./(q*p); %eq33
eq43=1-q.*A.*p*cos(2*pi*q) +q.*B.*p*sin(2*pi*q) - q^2.*p.*cos(phi)/(q^2-1);
eq44=q^2.*A*sin(2*pi*q) + q^2.*B*cos(2*pi*q) + q^2.*sin(phi)/(q^2-1);

% now get common roots
r1 = roots(eq43, eq44, 'ms'); % marching squares
r2 = roots(eq43, eq44); % figure out the method
r3 = roots(eq43,eq44, 'resultant'); % resulant
r4 = roots([eq43 ; eq44], 'ms'); % marching squares
r5 = roots([eq43 ; eq44]); % figure out the method
r6 = roots([eq43 ;eq44], 'resultant'); % resulant

pass(1) = norm( r1 - r2 ) < tol; 
pass(2) = norm( r1 - r3 ) < tol; 
pass(3) = norm( r1 - r4 ) < tol; 
pass(4) = norm( r1 - r5 ) < tol; 
pass(5) = norm( r1 - r6 ) < tol; 

end