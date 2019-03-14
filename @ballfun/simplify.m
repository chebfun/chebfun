function f = simplify( f, pref )
%SIMPLIFY  BALLFUN simplification
%
% F = SIMPLIFY( F ) returns a ballfun object simplified to have a
% compressed internal dimensions of coefficient tensor.
%
% This function is for internal use only.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( f )
    return
end

cfs = f.coeffs;
vals = ballfun.coeffs2vals( cfs );

vscl = max(1, max( abs( vals(:) ) ));

r_cfs = max(max( abs(cfs), [], 2), [], 3);
l_cfs = max(max( abs(cfs), [], 1), [], 3);
l_cfs = l_cfs(:);
t_cfs = max(max( abs(cfs), [], 1), [], 2);
t_cfs = t_cfs(:);

rTech = chebtech2.make( {'',r_cfs} );
lTech = trigtech.make( {'',l_cfs} );
tTech = trigtech.make( {'',t_cfs} );

rvals = rTech.coeffs2vals(rTech.coeffs);
rdata.vscale = vscl;
rdata.hscale = 1;
lvals = lTech.coeffs2vals(lTech.coeffs);
ldata.vscale = vscl;
ldata.hscale = 1;
tvals = tTech.coeffs2vals(tTech.coeffs);
tdata.vscale = vscl;
tdata.hscale = 1;

% Check happiness along each slice:
[resolved_r, cutoff_r] = happinessCheck(rTech, [], rvals, rdata);
[resolved_l, cutoff_l] = happinessCheck(lTech, [], lvals, ldata);
[resolved_t, cutoff_t] = happinessCheck(tTech, [], tvals, tdata);

% Simplify: 
if ( resolved_r )
    cfs = cfs(1:cutoff_r, :, :); 
end
n = size(cfs,2); 
mid = floor(n/2)+1;
if ( resolved_l )
    cfs = cfs(:, mid-floor(cutoff_l/2):mid+cutoff_l-floor(cutoff_l/2)-1, :); 
end
p = size(cfs,3); 
mid = floor(p/2)+1;
if ( resolved_t )
    cfs = cfs(:, :, mid-floor(cutoff_t/2):mid+cutoff_t-floor(cutoff_t/2)-1); 
end

f.coeffs = cfs;
end