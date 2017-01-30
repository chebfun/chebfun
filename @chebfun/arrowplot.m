function varargout = arrowplot(f,varargin)
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
%   A = []; for k = 1:3, A = [A exp((-.1*k+1i)*t)]; end, arrowplot(A,'color','k')
if ~isreal(f)
    if nargin < 2
        opts = {};
    end
else
    g = varargin{1};
    varargin(1) = [];
    f = f + 1i*g;
    if nargin < 3
        opts = {};
    end
end
fp = diff(f);
fdom = domain(f);
fend = feval(f, fdom(end));
fpend = feval(fp, fdom(end));
fpend = .001*fpend/norm(fpend,2);
pp = plot(f,varargin{:});
for aCounter=1:length(fend)
    h(aCounter)=annotation('arrow');
    set(h(aCounter),'parent', gca, ...
        'position', [real(fend(aCounter)) imag(fend(aCounter)) ...
            real(fpend(aCounter)) imag(fpend(aCounter))], ...
        'HeadLength', 10, 'HeadWidth', 10, 'HeadStyle', 'vback2', ...
        'Color', pp(aCounter).Color, varargin{:});
end

if ( nargout == 1)
    varargout{1} = pp;
end

end
