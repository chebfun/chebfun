function displayIVPinfo(u, isIVP, varargin)
%DISPLAYINFOIVP    Utility routine for displaying IVP solving
%
% Calling sequence:
%   DISPLAYINFO(VARARGIN)
% where
%   U:          The solution computed.
%   ISIVP:      Equal to 1 if we've just solved an IVP, 0 otherwise (FVP).
%   VARARGIN:   Various useful information passed from the Newton solver
%               to this method.
%
% See also: chebop/displayInfo.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developers note: This method was introduced to allow the IVP solver to use
% similar syntax as the BVP solver for displaying useful information about the
% solution process.
%
% At the moment, we only display information at the end of solving IVPs, rather
% than during the Newton iteration which is the case for BVPs. Thus, this method
% only gets called at the end of the solveivp method, which is why we don't need
% the 'mode' argument like we do in the BVP case.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Are we dealing with a CHEBFUN or a CHEBMATRIX?
if (isa(u, 'chebfun'))
    len = length(u);
else
    len = max(cellfun(@length, u.blocks));
end

% Did we just solve an IVP or an FVP?
if ( isIVP == 1 )
    fprintf('Initial value problem detected.\n');
else
    fprintf('Final value problem detected.\n');
end

% Show information about the number of pieces and the total length of the
% solution:
fprintf('Number of pieces of the solution: %i.\n', length(u.domain) - 1);
fprintf('Total length of solution: %i.\n', len);
   
end

