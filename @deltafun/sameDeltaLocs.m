function varargout = sameDeltaLocs(f, g)
%SAMEDELTALOCS  Finds common locations of delta functions in F and G.
%   idx1 = SAMEDELTALOCS(F, G) returns a logical vector of of same length as
%   that of F.DELTALOCS indicating which locations are shared with G.
%   [idx1, idx2] = SAMEDELTALOCS(F, G) does the same but also returns the
%   vector idx2 indicating locations of G.
%
% See also SIMPLIFY, CLEANROWS, CLEANCOLUMNS

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deat with the empty cases:
if ( isempty(f) || isempty(g) )
    if ( nargout == 2 )
        varargout{1} = [];
        varargout{2} = [];
    else
        varargout{1} = [];        
    end
    return;
end

% Initialize variables:
loc1 = f.deltaLoc;
idx1 = zeros(size(loc1));
loc2 = g.deltaLoc;
idx2 = zeros(size(loc2));

% Get preferences:
pref = chebfunpref();
tol = pref.deltaPrefs.proximityTol;

% Loop and compare locations:
for i = 1:length(loc1)
    for j = 1:length(loc2)
        % If both entries are equal to zero:
        p = ( loc1(i) == 0 && loc2(j) == 0 );
        % Or if they are very close to each other:
        maxLoc = max(abs([loc1(i), loc2(j)]));
        p = p | ( abs((loc1(i)-loc2(j)))/maxLoc < tol );    
        if ( p )
            % Mark the indices:
            idx1(i) = 1;
            idx2(j) = 1;
        end
    end
end

% Convert to logical indices:
idx1 = logical(idx1);
idx2 = logical(idx2);

% Output based on the number of outputs:
if ( nargout == 2 )
    varargout{1} = idx1;
    varargout{2} = idx2;
else
    varargout{1} = idx1;
end

end