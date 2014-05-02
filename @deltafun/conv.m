function h = conv(f, g)
%CONV   Convolution of DELTAFUN objects.
%   H = CONV(F, G) produces the convolution of DELTAFUN objects F and G:
%                     - 
%                    /
%           H(x) =   |    F(t) G(x-t) dt,  x in [a + c, b + d]
%                    /
%                   -

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Return empty for an empty input:
if ( isempty(f) || isempty(g) )
    h = deltafun();
    return
end

if ( ~isa(f, 'deltafun') )
    % Then g must be a deltafun
    h = conv(g, f);
    return
end
  
if ( ~isa(f, 'deltafun') )    
    deltaMagF = [];
    deltaLocF = [];
    funF = f;
else
    f = simplifyDeltas(f);
    deltaMagF = f.deltaMag;
    deltaLocF = f.deltaLoc; 
    funF = f.funPart;
end

if ( ~isa(g, 'deltafun') )    
    deltaMagG = [];
    deltaLocG = [];
    funG = g;
else
    g = simplifyDeltas(g);
    deltaMagG = g.deltaMag;
    deltaLocG = g.deltaLoc; 
    funG = g.funPart;
end


% Extract the domains of f and g:
domF = funF.domain;
domG = funG.domain;
a = domF(1);
b = domF(end);

c = domG(1);
d = domG(end);

% Get the threshold for deltafunctions:
pref = chebpref();
deltaTol = pref.deltaPrefs.deltaTol;

% Compute the convolution of funParts and append it to the output cell:
h = conv(funF, funG);


%% Get all the deltafunction contributions
% (f + df) * (g + dg) = f*g + dg * (f + df) + df * (g + dg) - df * dg
%                     = f*g + dg * (f + df/2) + df * (g + dg/2)

% Contributions due to deltafunctions in F:
% df * (g + dg/2):
if( ~isempty(deltaLocF))
    [m, n] = size(deltaMagF);
    % loop through the delta function matrix
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
                h = [h, {hij}]; %#ok<AGROW>
            end
        end
    end
end

% Contributions due to delta functions in G:
% dg * (f + df/2);
if ( ~isempty(deltaLocG) )
    [m, n] = size(deltaMagG);
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
                h = [h, {hij}]; %#ok<AGROW>                    
            end
        end
    end
end

end