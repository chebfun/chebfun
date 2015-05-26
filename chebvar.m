function chebvar(varargin)
%CHEBVAR   Short-cut for constructing CHEBFUN variables.
%   CHEBVAR arg1 arg2 ...
%   is short-hand notation for creating symbolic variables
%          arg1 = chebfun('arg1');
%          arg2 = chebfun('arg2'); ...
%   The outputs are created in the current workspace.
%
%   CHEBVAR arg1 arg2 ... DOM constructs the CHEBFUN objects on the domain DOM,
%   i.e.,  arg1 = chebfun('arg1', DOM);
%          arg2 = chebfun('arg2', DOM); ...
%
%   In both cases, the CHEBFUN is created according to the currently stored
%   default preferences.
%
%   Example:
%     chebvar x
%     f = sin(x)
%
% See also CHEBFUN, CHEBFUNPREF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Trivial case:
if ( nargin == 0 )
    return
end
dom = [];

% Locate valid variable names:
isVar = cellfun(@isvarname, varargin);

% The last entry may be a domain:
if ( ~all(isVar) )
    dom = str2num(varargin{end}); %#ok<ST2NM>
    if ( ~isempty(dom) )
        varargin(end) = [];
        isVar(end) = [];
    end
end

% Check validity of variable name:
if ( ~all(isVar) )
    error('CHEBFUN:chebvar:badName', 'Not a valid variable name.');
end

% Acquire some preferences:
pref = chebfunpref();
if ( isempty(dom) )
    dom = pref.domain;
end
    
% Loop over each of the inputs:
for k = numel(varargin):-1:1
    op = str2op(varargin{k});
    f = chebfun(op, dom, pref);
    assignin('caller', varargin{k}, f);
end

end

function op = str2op(op)
% This is here as it's a clean function with no other variables hanging around
% in the scope.
depVar = symvar(op);
if ( numel(depVar) ~= 1 )
    error('CHEBFUN:chebvar:indepVars', ...
        'Incorrect number of independent variables in string input.');
end
op = eval(['@(' depVar{:} ')', op]);
end
