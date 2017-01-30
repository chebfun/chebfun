function arrowplotNew(u,opts)
%ARROWPLOT   Chebfun plot with arrowhead at end
%   ARROWPLOT(F,G), where F and G are CHEBFUNs with the same
%   domain, plots the curve (F,G) in the plane with an arrowhead.
%
%   ARROWPLOT(F), where F is a complex CHEBFUN, plots the curve
%   (real(F),imag(F)) in the plane with an arrowhead.
%
%   Arguments can also be quasimatrices, in which case several
%   curves are plotted.
%
%   Plotting options can be passed in the usual fashion.
%
% Examples:
%
%   t = chebfun('t',[0,6]);
%   f = sin(t); g = cos(t), arrowplot(f,g);
%
%   h = exp((-.2+3i)*t); arrowplot(h,'color','r')
%
%   A = []; for k = 1:3, A = [A exp(-.1*k+1i)*t]; arrowplot(A) 
if nargin < 2
    opts = {};
end
f = u + 1i*diff(u);
fp = diff(f);
fend = f(end);
fpend = fp(end);
fpend = .2*fpend/norm(fpend,2);
plot(f,opts{:});
h=annotation('arrow');
set(h,'parent', gca, ...
    'position', [real(fend) imag(fend) real(fpend) imag(fpend)], ...
    'HeadLength', 10, 'HeadWidth', 10, 'HeadStyle', 'vback2', ...
    opts{:});

end
