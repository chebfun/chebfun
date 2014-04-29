function out = domainCheck(f, g)
%DOMAINCHECK   True if the domains of two CHEBFUN2 objects are the same.
%   DOMAINCHECK(F, G) returns TRUE if the domains of the two
%   CHEBFUN2 objects F and G coincide up to a tolerance depending on their
%   horizontal scales or if both F and G are empty CHEBFUN objects.
%
% See also CHEBFUN/DOMAINCHECK.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The CHEBFUN2 class uses this function internally to compare the domains of
% CHEBFUN2 objects before attempting to perform operations on multiple CHEBFUN
% objects that require the CHEBFUN2 objects to reside on the same interval.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Empty check: 
if ( isempty( f ) && isempty( g ) ) 
    out = true; 
    return
elseif ( xor(isempty( f ), isempty( g ) ) )
    out = false; 
    return
end

% Extract the columns and rows: 
fcols = f.cols; 
frows = f.rows; 
gcols = g.cols; 
grows = g.rows; 

% Call 1D domain check: 
out = domainCheck(fcols, gcols) & domainCheck(frows, grows); 

end
