function arrowplotNew(u,opts)  % plot chebfun with arrow on end
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
