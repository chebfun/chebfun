function out = vscale(F, supString)
%VSCALE   Vertical scale of a CHEBFUN object.
%   VSCALE(F) returns an estimate of the maximum absolute value of F. VSCALE
%   always returns a scalar, even when F is an array-valued CHEBFUN or a
%   quasimatrix. Vertical scales of each of the piecewise components and columns
%   of F are given by get(F, 'vscale-local'); Values of F at its break points
%   are ignored by VSCALE(F).
% 
%   VSCALE(F, 'ess-sup') is the same as VSCALE(F)
% 
%   VSCALE(F, 'sup') also takes into account the point values of the object 
%   F at its break points while computing its VSCALE. This is the vscale
%   returned in CHEBFUN/DISPLAY.
%
% See also MAX, MINANADMAX.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(F) )
    out = [];
    return
end

if ( nargin < 2 )
    % By default, we ignore values at the break points while computing the
    % vsclase:
    ignoreBreaks = 1;    
else
    switch supString
        case 'ess-sup'          % This is the same as default.
            ignoreBreaks = 1;
        case 'sup'
            ignoreBreaks = 0;   % Don't ignore break point values.
        otherwise
            error( 'CHEBFUN:CHEBFUN:vscale:unknownFlag', 'unknown flag passed.' )
    end
end
            
out = 0;
for k = 1:numel(F)
    % Get the local vscales:
    v = get(F(k), 'vscale-local');

    % Compute the maximum:
    out = max(out, max(v(:)));
    if ( ~ignoreBreaks )
        out = max(out, max(abs(F(k).pointValues(:))));
    end
end

end
