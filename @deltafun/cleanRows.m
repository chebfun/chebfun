function A = cleanRows(A)
% Remove trivial rows:

tol = deltafun.pref.deltafun.deltaTol;

while( max(abs(A(end, :))) < tol )
    A(end, :) = [];
end