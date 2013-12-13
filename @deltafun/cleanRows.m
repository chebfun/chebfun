function A = cleanRows(A)
% Remove trivial rows:

tol = deltafun.pref.deltafun.deltaTol;

if ( isempty(A) )
    return
end

while( max(abs(A(end, :))) < tol )
    A(end, :) = [];
end