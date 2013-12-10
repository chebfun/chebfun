function s = plus(f, g)
%+   Addition of DELTAFUN objects.
%   F + G adds F and G, where F and G may be DELTAFUN objects or scalars.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%
s = deltafun;
% If one of the arguments is empty:
if ( isempty(f) || isempty(g) )
    % return with the empty DELTAFUN:
    return;
end

% One of the arguments i.e. f or g is necessarily a DELTAFUN object. Otherwise, 
% this overloaded plus would not have been called.

% If g is a double, f is a deltafun:
if ( isa(g, 'double') )
    s = f;
    s.funPart = s.funPart + g;
end

% Recursive call, by making the deltafun g as the first argument:
if ( isa(f, 'double') )
    s = g + f;
    return
end

if ( isa(f, 'deltafun') && isa(g, 'deltafun') )
    %[TODO]: This should be based on tolerances?
    if ( f.domain ~= g.domain ) 
        error( 'CHEBFUN:DELTAFUN:plus', 'f and g must have the same domain' );
    else
        s.domain = f.domain;
    end
    
    %[TODO]: what shoudl be here?
    s.isTransposed = 0;
    
    s.funPart = f.funPart + g.funPart;
    
    % Add the delta functions:
    
    va = f.delta.location;
    vb = g.delta.location;
    
    A = f.delta.magnitude;
    B = g.delta.magnitude;
    
    [C, vc] = sortAndMerge(A, va, B, vb);

    s.delta.location = vc;
    s.delta.magnitude = C;
    
end

%% %%%%%%%%%%%%% DELTAFUN + CHEBFUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This should be removed, since CHEBFUNs should be casted to DELTAFUNs first and
% then added
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recursive call: If g is a chebfun, upgrade it to a deltafun and call plus
% again.
if ( isa(f, 'deltafun') && isa(g, 'chebfun') )
    %[TODO]: This should be based on tolerances?
    if ( f.domain ~= g.domain ) 
        error( 'CHEBFUN:DELTAFUN:plus', 'f and g must have the same domain' );
    end
    t = deltafun.zeroDeltaFun(g.domain);
    t.funPart = g;
    s = f + t;
    return
end

%%
% If the result is a smooth object, return a smooth object.
if ( issmooth(s) )
    % return a smooth object if the operation has removed all 
    % delta functions.
    s = s.funPart;
end

end
%%
function [C, vc] = sortAndMerge(A, va, B, vb)
% SortAndMerge sorts and merges matrices A and B based the values in va and vb.
% The output is a matrix C which contains columns of A and B merged according to
% the sorted values in the combined set {va, vb}. Duplication in {va, vb}
% results in the corresponding columns being added = merged into a single
% column. It is assumed that va and vb have unique elements individually.

% Get the tolerance
tol = deltafun.pref.deltafun.proximityTol;

% Make the number of rows in A and B same by appending zero rows.
ra = size(A, 1); 
rb = size(B, 1);
if( ra >= rb )
    B = [ B; zeros(ra-rb, size(B, 2)) ];
else
    A = [ A; zeros(rb-ra, size(A, 2)) ];
end

% Make sure the input is sorted:
[va, idx] = sort(va);
A = A(:, idx);

[vb, idx] = sort(vb);
B = B(:, idx);


% Initialize variables:
m = 1; n = 1; i = 1;
lenA = length(va);
lenB = length(vb);

vc = zeros(1, lenA + lenB);
C  = zeros(size(A,1), length(vc)); % or zeros(size(B, 1))

while ( m <= lenA || n <= lenB )
    % If one of the arrays is exausted, copy the second one into output, no
    % duplication is assumed in va or vb, so:
    if ( m > lenA )
        while ( n <= lenB )
            vc(i) = vb(n);
            C(:, i) = B(:, n);
            n = n + 1;
            i = i + 1;
        end
        break;
    end
    
    if ( n > lenB )
        while ( m <= lenA )
            vc(i) = va(m);
            C(:, i) = A(:, n);
            m = m + 1;
            i = i + 1;
        end
        break;
    end
    
    % None of the arrays is exausted;
    
    % Duplication scenarios
    p = ( va(m) == 0 && vb(n) == 0 );
    vcMax = max(abs([va(m), vb(n)]));
    p = p | ( abs((va(m)-vb(n)))/vcMax < tol );
    % If there is a duplicate
    if ( p )
        C(:,i) = A(:, m) + B(:, n);
        vc(i) = va(m);
        m = m + 1;
        n = n + 1;
    elseif ( va(m) < vb(n) )
        vc(i) = va(m);
        C(:, i) = A(:, m);
        m = m + 1;
    else
        vc(i) = vb(n);
        C(:, i) = B(:, n);
        n = n + 1;
    end
    i = i + 1;
end
% compensate for the last spurious addition in i and truncate the matrices to
% appropriate lengths.
C = C(:, 1:i-1);
vc = vc(1:i-1);
end