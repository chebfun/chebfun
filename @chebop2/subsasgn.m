function varargout = subsasgn(f, index, varargin)
%SUBSASGN   SUBSASGN for CHEBOP2 objects. 
%
% Note: Designed so boundary conditions of an operator can be converted to 
% CHEBFUNs before a call to MLDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

idx = index(1).subs;
vin = varargin{:};

switch index(1).type
    
    case '.' 
        if ischar(idx)
            varargout = {set(f, idx, vin)};
        else
            error('CHEBFUN:CHEBOP2:subsasgn:UnexpectedType', ...
                  'Chebop2 properties are string arrays.')
        end   
end
end