function f = simplify(f)
%SIMPLIFY   Removes trivial rows of a DELTAFUN F.
%   SIMPLIFY(F) removes trivial rows from the magnitude matrix of the 
%   DELTAFUN F based on the tolerance. 
%
% See also SUM, CUMSUM.
deltaMag = f.delta.magnitude;
while( max(abs(deltaMag(end,:))) < deltafun.pref.deltafun.deltaTol )
    deltaMag = deltaMag(1:end-1,:);
end 
f.delta.magnitude = deltaMag;
