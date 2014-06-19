function pass = test_chap21(pref)

% Here we test that everything in chapter 21 of ATAP runs OK.

warnState = warning('off', 'CHEBFUN:CHEBOP:feval:deprecated');

x = chebfun('x'); p = sin(x); length(p);
pp = diff(p); x14 = chebpts(14); pp14 = pp(x14);
D = chebop(@(u) diff(u)); D14 = D(14);
pass(1) = norm(pp14-D14*p(x14)) < 1e-10;
D(3);
D(5);
D2 = chebop(@(u) diff(u,2)); D2(5);
norm(D2(5)-(D(5))^2);
p3 = chebfun([0 0 0 1 0]'); FS = 'fontsize';
% clf, plot(p3,'.-'), title('Lagrange polynomial l_3',FS,9)
p3pp = diff(p3,2); x5 = chebpts(5); p3pp(x5);
pass(2) = all(size(D(33)) == 33);
f = sin(7*x).*exp(x).*tan(x); norm(diff(f)-D*f);
L = chebop(@(u) diff(u,2) + diff(u) + 100*u);
L(5);
f = exp(x); Lf = L*f;
Lfexact = 102.*exp(x); 
pass(3) = norm(Lf-Lfexact) < 1e-10;
L.bc = 'dirichlet'; feval(L,5,'oldschool');
x5 = chebpts(5); x5([1 end]) = 0;
u5 = feval(L,5,'oldschool')\x5;
% plot(chebfun(u5),'.-')
% title('Spectral solution to (21.3) on 5-point grid',FS,9)
x12 = chebpts(12); x12([1 end]) = 0;
u12 = feval(L,12,'oldschool')\x12;
% plot(chebfun(u12),'.-')
% title('Spectral solution to (21.3) on 12-point grid',FS,9)
u = L\x; 
% plot(u,'.-')
% title(['Spectral solution to (21.3) on ' ...
%     'automatically determined grid'],FS,9)
length(u);
pass(4) = norm(L*u-x) < 1e-10;
L.bc = 'neumann'; feval(L,5,'oldschool');
u = L\x; 
% plot(u), ylim([-0.015 0.015])
% title('Solution to (21.3) except with Neumann BCs',FS,9)
L = chebop(@(x,u) diff(u,2)-x.*u,[-30,30]);
L.lbc = 1; L.rbc = 0; u = L\0; 
% plot(u), title('Solution to Airy equation (21.4)',FS,9)
N = chebop(0,6);
N.op = @(theta) diff(theta,2) + sin(theta);
N.lbc = -pi/2; N.rbc = pi/2; theta = N\0;
% plot(-cos(theta)), grid on, ylim([-1 1])
% title('Nonlinear pendulum (21.5)',FS,9)
% xlabel('t',FS,10), ylabel('height -cos(\theta)')
N.lbc = -pi/2; N.rbc = 5*pi/2; theta = N\0;
% plot(-cos(theta)), grid on, ylim([-1 1])
% title('Nonlinear pendulum (21.5), another solution',FS,9)
% xlabel('t',FS,10), ylabel('height -cos(\theta)')

warning(warnState);

end
