% Tests for chebfun plotting functions.
function pass = test_plot_xylim(pref)

tol = 1e-4;
% Create a figure, and make it invisible. Need to do this a number of time
% throughout the test.
hfig = figure;
set(hfig,'visible','off')
%% Finite functions on unbounded domains

dom1 = [0 pi];
pass1 = [];
x = chebfun(@(x) x, dom1);
plot(x)
pass1(length(pass1) + 1) = ( norm(dom1 - get(gca,'xlim')) < tol);

plot(0.62*sin(x))
yl = get(gca, 'ylim');
pass1(length(pass1) + 1) = ( norm(0.62 - yl(2)) > 0.05);

plot(0.62*[sin(x) 0*x -sin(x)])
yl = get(gca, 'ylim');
pass1(length(pass1) + 1) = ( norm(0.62*[-1 1] - yl) > 0.05);
pass1(length(pass1) + 1) = strcmp(get(gca,'ylimmode'), 'auto');

plot(sin(x))
hold on
plot(-sin(x))
% yl = get(gca,'ylim');
% pass1(length(pass1) + 1) = ( norm(yl - [-1 1]) < tol );
pass1(length(pass1) + 1) = strcmp(get(gca,'ylimmode'), 'auto');
hold off

plot(.62*sin(x))
hold on
plot(-.62*sin(x))
yl = get(gca,'ylim');
pass1(length(pass1) + 1) = (norm(0.62*[-1 1] - yl) > 0.05 );
pass1(length(pass1) + 1) = strcmp(get(gca,'ylimmode'), 'auto');
hold off

%% Unbounded functions

dom2 = [-inf 0];
pass2 = [];
f = chebfun(@(x) exp(x), dom2);
plot(f)
pass2(length(pass2) + 1) =  strcmp(get(gca,'xlimmode'), 'manual');
pass2(length(pass2) + 1) =  strcmp(get(gca,'ylimmode'), 'auto');
xl = get(gca,'xlim');
pass2(length(pass2) + 1)  = ( norm(xl - [-10 0]) < tol );

g = chebfun(@(x) -0.62*exp(x), dom2);
plot(g)
pass2(length(pass2) + 1) =  strcmp(get(gca,'ylimmode'), 'auto');
yl = get(gca,'ylim');
pass2(length(pass2) + 1)  = ( norm(yl(1) - 0.62) > 0.05 );

dom3 = [-20 20];
h = chebfun(@(x) cos(x), dom3);
hold on
plot(h,'r')
xl = get(gca,'xlim');
pass2(length(pass2) + 1) = ( norm(xl - dom3) < tol );
hold off

plot(h,'g')
hold on
plot(g)
xl = get(gca,'xlim');
pass2(length(pass2) + 1)  = ( norm(xl - dom3) < tol );
hold off

dom4 = [-2 2];
h = chebfun(@(x) cos(x), dom4);
plot(g)
hold on
plot(h,'r')
xl = get(gca,'xlim');
pass2(length(pass2) + 1) = ( norm(xl - [-10 2]) < tol );
hold off

plot(h,'g')
hold on
plot(g)
xl = get(gca,'xlim');
pass2(length(pass2) + 1)  = ( norm(xl - [-10 2]) < tol );
hold off
%% Functions that blow-up
pass3 = [];
f1 = 11.3*sin(x);
f2 = 1./x;

plot(f2)
pass3(length(pass3) + 1) =  strcmp(get(gca,'xlimmode'), 'manual');
pass3(length(pass3) + 1) =  strcmp(get(gca,'ylimmode'), 'manual');
xl = get(gca, 'xlim');
yl = get(gca, 'ylim');
% Check that we obtain reasonable ylimits
pass3(length(pass3) + 1)  = ( norm(yl(1) - 1/xl(2)) < tol && ( yl(2) < 10 ) );
hold on
plot(f1,'r')
% Check that we obtain reasonable ylimits
yl = get(gca, 'ylim');
pass3(length(pass3) + 1)  = ( norm(yl(1)) < tol && ( yl(2) > 10 ) );
hold off

% Do the plotting in reverse order
plot(f1,'g')
ylOld = get(gca, 'ylim');
hold on
plot(f2)
% Check that we obtain reasonable ylimits
ylNew = get(gca, 'ylim');
pass3(length(pass3) + 1) = ( (norm(ylOld(1) - ylNew(1)) < tol) && ...
    (ylNew(2) > 10) );
hold off

% Check x-lim and y-lim symmetry of x:
f = chebfun(@(x) x, [-inf inf]);
plot(f)
xl = get(gca, 'xlim');
yl = get(gca, 'ylim');
pass3(length(pass3) + 1) = (abs(sum(xl)) < 1e-10) && (abs(sum(yl)) < 1e-10);

%%
pass = [pass1, pass2, pass3];

close(hfig);
