function pass = test_chebmatrix_spy(pref)

pass(1) = true;

return

%%

close all
dom = [-1 0 1];
I = operatorBlock.eye(dom);
D = operatorBlock.diff(dom);
x = chebfun('x', dom);
c = sin(x.^2);
C = operatorBlock.mult(c);   
E = functionalBlock.eval(dom);
El = E(dom(1));
L = [ D^2, -I, sin(x); C, D, chebfun(0,dom); functionalBlock.zero(dom), El, 4 ];
spy(L, 'disc', 'colloc2')
open test_chebmatrix_spy1.fig
shg

%%

close all
spy(L, 'disc', 'ultraS')
open test_chebmatrix_spy2.fig
shg

end

