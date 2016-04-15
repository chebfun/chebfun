function out = domainCheck(f, g)
%DOMAINCHECK   True if the domains of two CHEBFUN3 objects are the same.
%   DOMAINCHECK(F, G) returns TRUE if the domains of the two
%   CHEBFUN3 objects F and G coincide up to a tolerance depending on their
%   horizontal scales or if both F and G are empty CHEBFUN objects.
%
% See also CHEBFUN/DOMAINCHECK.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The CHEBFUN3 class uses this function internally to compare the domains of
% CHEBFUN3 objects before attempting to perform operations on multiple
% CHEBFUN3 objects that require the CHEBFUN3 objects to reside on the same interval.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Empty check: 
if ( isempty(f) && isempty(g) )
    out = true;
    return
elseif ( xor(isempty(f), isempty(g) ) )
    out = false;
    return
end

% Extract the columns, rows and tubes: 
fcols = f.cols;
frows = f.rows;
ftubes = f.tubes;
gcols = g.cols;
grows = g.rows;
gtubes = g.tubes;

% Call 1D domain check: 
out = domainCheck(fcols, gcols) & domainCheck(frows, grows) & ...
    domainCheck(ftubes, gtubes);

end