function [ishappy, epslevel, cutoff] = looseCheck(f, values, pref)
%LOOSECHECK   Loose happiness check for Chebtech construction.
%   [ISHAPPY, EPSLEVEL, CUTOFF] = LOOSECHECK(F, PREF). This file is not yet
%   implemented.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Absolute value of coefficients:
ac = max(abs(f.coeffs(end:-1:1,:))./max(abs(values), [], 1), [], 2);

% Cumulative maximum:
cac = cummax(ac);

% Window it:
n = length(cac);
skip = round(min(5, n/4));
window = round(min(10, n/2));
nn = 5:skip:length(cac)-window;
q = nan(length(nn),1);
k = 0;
for m = nn
    k = k+1;
    q(k) = mean(cac(m:m+window));
end
q = q(end:-1:1);

if ( numel(q) < 5 )
    ishappy = false;
    epslevel = eps;
    cutoff = inf;
    return
end

lq = log10(q);
dq2 = diff(lq, 2);
dq2(dq2 < 0) = 0;

% The magix formula:
r = abs(lq(3:end).^5).*(dq2);
[ignored, idx] = max(r);
cutoff = nn(idx+2);
qq = q(idx:end);
epslevel = max(mean(qq)/5, min(qq));
tol = max(n*(pref.eps).^(2/3)/epslevel, 100*eps);
grad = abs(mean(diff(q(idx+1:end))));
ishappy = r(idx) > 1 && grad < tol;

% % Plotting for testing:
% f.epslevel = epslevel;
% chebpolyplot(f), hold on
% semilogy(nn,abs(q),'r'), shg, hold off
% figure(2)
% plot(nn(3:end),r,'-m');
% figure(1)
% pause

if ( isempty(ishappy) )
    ishappy = false;
    epslevel = eps;
end

end

function y = cummax(a)
n = size(a,1);
y = zeros(size(a)); 
a = a(end:-1:1);
y(1) = a(1); 
for i=2:n 
   y(i) = max( y(i-1), a(i) ); 
end 
end
