clc
clear classes

%% Building blocks
dom = [-2 2];
I = linop.eye(dom);
D = linop.diff(dom);
Z = linop.zeros(dom);
x = chebfun('x', dom);
c = sin(x.^2);
C = linop.diag(c);   
E = linop.eval(dom);
El = E(dom(1));
Er = E(dom(end));


%% Solve a linear system 
L = linop([ D^2, -I, x; C, D, 0*x; linop.zero(dom), El, 4 ]);
f = [abs(x-1); 0*x; 1 ];
B1 = [El, -Er, 0];
B2 = [linop.sum(dom), El, 0];
B3 = [Er*D, linop.zero(dom), 0];
B4 = [Er,-El,2];
L = addbc(L,B1,0);
L = addbc(L,B2,1);
L = addbc(L,B3,0);
%L = addbc(L,B4,0);

%%
w = L\f;

%%
plot(w{1},'b'); hold on
plot(w{2},'r'); hold off, shg
w{3}

%%
% check the ODEs
zero1 = norm( diff(w{1},2)-w{2}+x*w{3} - f{1} )
zero2 = norm( c.*w{1} + diff(w{2}) + 0 - f{2} )
zero3 = abs( 0 + feval(w{2},dom(1)) + 4*w{3} - f{3} )

%%
% check the BCs
v = w{2};  u = w{1};
zero4 = abs( u(-2)-v(2) )
zero5 = abs( sum(u)+v(-2) - 1)
zero6 = abs( feval(diff(u),dom(end)) )
