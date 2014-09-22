classdef nonPeriodicTech < smoothfun % (Abstract) 
%NONPERIODICTECH   Approximate smooth functions on [-1,1].
%   Abstract (interface) class for approximating functions on the 
%   interval [-1,1], with a basis of nonperiodic functions.
%
% See also SMOOTHFUN, CHEBTECH.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NONPERIODICTECH Class Description:
%
% The NONPERIODICTECH class is an abstract class for representations of 
% functions on the interval [-1,1], with a basis of nonperiodic functions.
%
% The current instance of PERIODICTECH is CHEBTECH.
%
% Class diagram: [<<SMOOTHFUN>>] <-- [<<NONPERIODICTECH>>] <-- [<<CHEBTECH>>]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS IMPLEMENTED IN THIS M-FILE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function out = isPeriodicTech(f)
        %ISPERIODIC    Test if the objtect is is constructed with a basis of
        %periodic functions. 
        %    Returns 0 for NONPERIODICTECH.
            out = 0;
        end
        
    end
    
end
