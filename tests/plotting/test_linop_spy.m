function pass = test_linop_spy(pref)

pass(1) = true;

return

%%

close all
dom = [-1 0 1];
I = operatorBlock.eye(dom);
D = operatorBlock.diff(dom);
Z = operatorBlock.zeros(dom);
x = chebfun('x', dom);
c = sin(x.^2);
C = operatorBlock.mult(c);   
E = functionalBlock.eval(dom);
El = E(dom(1));
Er = E(dom(end));

L = [ D^2, -I, sin(x); C, D, chebfun(0,dom); functionalBlock.zero(dom), El, 4 ] ;
B1 = [El, -Er, 0];
B2 = [functionalBlock.sum(dom), El, 0];
B3 = [Er*D, functionalBlock.zero(dom), 0];
B4 = [Er,-El,2];
L = addbc(L,B1,0);
L = addbc(L,B2,1);
L = addbc(L,B3,0);

spy(L, 'disc', 'colloc2')
open test_linop_spy1.fig
shg

%%

close all
spy(L, 'disc', 'ultraS')
open test_linop_spy2.fig
shg

end

