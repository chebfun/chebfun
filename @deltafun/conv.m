function h = conv(f, g)
%CONV   Convolution of DELTAFUN objects.
%   H = CONV(F, G) produces the convolution of DELTAFUN objects F and G:
%                     - 
%                    /
%           H(x) =   |    F(t) G(x-t) dt,  x in [a + c, b + d]
%                    /
%   The output H is a cell array of DELTAFUNs or CLASSICFUNS. If the result of
%   the convolution of F and G is zero, an empty cell is returned. The cell
%   array returned is used by higher level convolutions.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Return empty for an empty input:
if ( isempty(f) || isempty(g) )
    h = deltafun();
    return
end

if ( ~isa(f, 'deltafun') )
    deltaMagF = [];
    deltaLocF = [];
    funF = f;
else
    f = simplifyDeltas(f);
    if ( isa(f, 'deltafun') )
        deltaMagF = f.deltaMag;
        deltaLocF = f.deltaLoc; 
        funF = f.funPart;
    else
        deltaMagF = [];
        deltaLocF = [];
        funF = f;
    end
end

if ( ~isa(g, 'deltafun') )
    deltaMagG = [];
    deltaLocG = [];
    funG = g;
else
    g = simplifyDeltas(g);
    if ( isa(g, 'deltafun') )
        deltaMagG = g.deltaMag;
        deltaLocG = g.deltaLoc; 
        funG = g.funPart;
    else
        deltaMagG = [];
        deltaLocG = [];
        funG = g;
    end
end

% Extract the domains of f and g:
domF = funF.domain;
domG = funG.domain;
a = domF(1);
b = domF(end);
c = domG(1);
d = domG(end);

% Get the threshold for deltafunctions:
pref = chebfunpref();
deltaTol = pref.deltaPrefs.deltaTol;

% Compute the convolution of funParts and append it to the output cell:
h = conv(funF, funG);

% Remove zero funs for simplicity:
zeroIndices = cellfun( @(hk) iszero(hk), h);
h(zeroIndices) = [];

%% Get all the deltafunction contributions
% (f + df) * (g + dg) = f*g + dg * (f + df) + df * (g + dg) - df * dg
%                     = f*g + dg * (f + df/2) + df * (g + dg/2)

% Contributions due to deltafunctions in F:
% df * (g + dg/2):
[m, n] = size(deltaMagF);
% Loop through the delta function matrix:
for i = 1:m
    for j = 1:n
        if ( abs(deltaMagF(i, j)) > deltaTol )
            % Take appropriate derivative, scale and shift the function:
            hij = deltaMagF(i, j) * changeMap(diff(g,i-1), deltaLocF(j) + [c d]);
            % The half below is to make sure that delta-delta interaction
            % is not counted twice:                
            if ( isa(hij, 'deltafun') )
                hij.deltaMag = hij.deltaMag/2;
            end

            % If the result is non-zero, append it:
            if ( ~iszero(hij) )
                h = [h, {hij}]; %#ok<AGROW>
            end
        end
    end
end

% Contributions due to delta functions in G:
% dg * (f + df/2);
[m, n] = size(deltaMagG);
% Loop through the delta function matrix:    
for i = 1:m
    for j = 1:n
        if ( abs(deltaMagG(i, j)) > deltaTol )
            % Take appropriate derivative, scale and shift the function
            hij = deltaMagG(i, j) * changeMap(diff(f,i-1), deltaLocG(j) + [a b]);
            % The half below is to make sure that delta-delta interaction
            % is not counted twice:
            if ( isa(hij, 'deltafun') )
                hij.deltaMag = hij.deltaMag/2;
            end

            % If the result is non-zero, append it:
            if ( ~iszero(hij) )
                h = [h, {hij}]; %#ok<AGROW>
            end                
        end
    end
end

end