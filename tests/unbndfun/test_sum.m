% function pass = test_sum(f, pref)
% 
% if ( nargin == 1 )
%     pref = fun.pref();
% end

f = unbndfun(@(x) exp(-x.^2), [-inf inf]);
sum(f) - sqrt(pi)

% f = unbndfun(@(x) x.^2.*exp(-x.^2), [-inf inf]);
% sum(f) - sqrt(pi)/2

% f = unbndfun(@(x) exp(-x.^2), [0 inf]);
% sum(f) - sqrt(pi)/2
% 
% f = unbndfun(@(x) exp(-x), [0 inf]);
% sum(f) - 1
% 
% f = unbndfun(@(x) x.*exp(-x), [0 inf]);
% sum(f) - 1

