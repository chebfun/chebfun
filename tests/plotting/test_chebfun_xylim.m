% Tests for chebfun plotting functions.
function pass = test_chebfun_ylim(pref)

tol = 1e-10;

%% Finite functions on unbounded domains
dom = [0 pi];
x = chebfun(@(x) x, dom);
plot(x)
pass(1) = ( norm(dom - get(gca,'xlim')) < tol);

plot(0.62*sin(x))
yl = get(gca, 'ylim');
pass(2) = ( norm(0.62 - yl(2)) > 0.05);

plot(0.62*[sin(x) 0*x -sin(x)])
yl = get(gca, 'ylim');
pass(3) = ( norm(0.62*[-1 1] - yl) > 0.05);
pass(4) = strcmp(get(gca,'ylimmode'), 'auto');

plot(sin(x))
hold on
plot(-sin(x))
yl = get(gca,'ylim');
pass(5) = ( norm(yl - [-1 1]) < tol );
pass(6) = strcmp(get(gca,'ylimmode'), 'auto');
hold off

%% Unbounded functions
dom2 = [-inf 0];
f = chebfun(@(x) exp(x), dom2);
plot(f)
pass(7) =  strcmp(get(gca,'xlimmode'), 'manual');
pass(8) =  strcmp(get(gca,'ylimmode'), 'auto');
xl = get(gca,'xlim');
yl = get(gca,'ylim');
pass(9)  = ( norm(xl - [-10 0]) < tol );
pass(10) = ( norm(yl - [0 1]) < tol );

g = chebfun(@(x) -0.62*exp(x), dom2);
plot(g)
pass(11) =  strcmp(get(gca,'ylimmode'), 'auto');
yl = get(gca,'ylim');
pass(12)  = ( norm(yl(1) - 0.62) > 0.05 );
shg

dom3 = [-20 20];
h = chebfun(@(x) cos(x), dom3);
hold on
plot(h,'r')
xl = get(gca,'xlim');
pass(13)  = ( norm(xl - dom3) < tol );
hold off

plot(h,'g')
hold on
plot(g)
xl = get(gca,'xlim');
pass(14)  = ( norm(xl - dom3) < tol );
hold off

dom4 = [-2 2];
h = chebfun(@(x) cos(x), dom4);
plot(g)
hold on
plot(h,'r')
xl = get(gca,'xlim');
pass(15)  = ( norm(xl - [-10 2]) < tol );
hold off

plot(h,'g')
hold on
plot(g)
xl = get(gca,'xlim');
pass(16)  = ( norm(xl - [-10 2]) < tol );
hold off

%% Functions that blow-up

f1 = 11.3*sin(x);
f2 = 1./x;

plot(f2)
pass(17) =  strcmp(get(gca,'xlimmode'), 'manual');
pass(18) =  strcmp(get(gca,'ylimmode'), 'manual');
xl = get(gca, 'xlim');
yl = get(gca, 'ylim');
pass(19)  = ( yl(1) < 1/xl(1) && ( yl(2) < 10 ) );
% pass(10) = ( norm(yl - [0 1]) < tol );

pass
