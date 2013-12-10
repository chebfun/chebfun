function prefs = pref(varargin)
%PREF   Preference settings for SINGFUN. [TODO:] SINGFUN will have no 
%   perefernces, instead the abstract superior class ONEFUN will have all the
%   preferences in future.
% 
%   SINGFUN.PREF(PREFNAME) returns the value corresponding to the preference
%   named in the string PREFNAME.
%
%   P = SINGFUN.PREF returns a structure P with a field P.SINGFUN which contains
%   the default SINGFUN preferences as fields/values. This structure may be
%   passed to the SINGFUN constructor.
%
%   P = SINGFUN.PREF(P) will check to see whether the input preference structure
%   P already has a SINGFUN field. If it does not, one is appended.
%
%   P = SINGFUN.PREF('PREFNAME1', VAL1, 'PREFNAME2', VAL2, ...) returns the same
%   structure as above, but with the default SINGFUN preferences 'PREFNAME1',
%   'PREFNAME2', replaced by the values in VAL1, VAL2, etc.
%
%   P = SINGFUN.PREF(P, 'PREFNAME1', VAL1, 'PREFNAME2', VAL2, ...) appends a
%   SINGFUN preference field to P if required, and modifies the SINGFUN
%   properties 'PREFNAME1' and 'PREFNAME2'.
%
%   Note that no checking of either the input PREFNAMEs or VALs takes place.
%
%   SINGFUN PREFERENCES (case sensitive)
%
%   [TODO] Copied from CHEBTECH, not done yet. Some of these might make sense
%   for SINGFUNs:
%
%     vsclae       -  What would they mean for a delta? 
%     hscale       -  //ditto//
%     deltaTol     -  Tolerance for delta functio magnitude.
%     [1e-12]      
%     maxDiffOrder -  Maximum Order of the derivative allowed for a deltafun
%       [20]
%                  

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

classname = 'deltafun';

% Has a preference structure been passed?
if ( numel(varargin) > 0 && isstruct(varargin{1}) ) % Yes
    prefs = varargin{1};
    varargin(1) = [];
else                                                % No
    prefs = struct();
end

% If it was, did it have a DELTAFUN field?
if ( isfield(prefs, classname) )  % It does, so either:
    if ( numel(varargin) == 0 )
        return                    % a) No props to change, return
    else
        p = prefs.(classname);    % b) Grab DELTAFUN prefs
    end
else
    % No DELTAFUN prefs found, so make some:
    p.deltaTol      = 1*1e-12;
    p.proximityTol  = 1*1e-12;
end
% p is now the preference substructure relating to the FUN class.

if ( isfield(prefs,'misc') ) 
    q = prefs.misc;
else
    q = struct();
end
% q is now the preference substructure for MISC preferences.

for names = fieldnames(q)'
    if isfield(p,names)
        field = names{:};
        p.(field) = q.(field);
        q = rmfield(q,field);
    end
end

% Two preference structures were passed. Copy matching fieldnames.
if ( numel(varargin) > 0 && isstruct(varargin{1}) ) % Yes
    r = varargin{1};
    varargin(1) = [];
    for names = fieldnames(r).'
        if isfield(p,names)
            field = names{:};
            p.(field) = r.(field);
        end
    end
end

% A single property has been queried, so return this property
if ( numel(varargin) == 1 )
    prefs = p.(varargin{1});
    return
end

% Property names have been passed, so alter/add DELTAFUN/MISC properties.
for k = 1:2:numel(varargin)
    if ( isfield(p, varargin{k}) )
        p.(varargin{k}) = varargin{k+1};
    else
        q.(varargin{k}) = varargin{k+1};
    end
end

% Append DELTAFUN preferences to the preference structure prefs for output.
prefs.(classname) = p;

% Append MISC preferences to the preference structure prefs for output.
prefs.misc = q;

end
