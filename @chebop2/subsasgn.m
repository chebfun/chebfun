function varargout = subsasgn(f,index,varargin)
% Subsasgn for a chebop2.  Designed so boundary conditions of an operator
% can be converted to chebfuns before a call to mldivide.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

idx = index(1).subs;
vin = varargin{:};

switch index(1).type
    case '.'
        if strcmp(idx,'bc')
            f.lbcshow = vin;
            f.rbcshow = vin;
            f.ubcshow = vin;
            f.dbcshow = vin;
            varargout = {set(f,'bc',vin)};
        elseif strcmp(idx,'lbc')
            f.lbcshow = vin;
            varargout = {set(f,'lbc',vin)};
        elseif strcmp(idx,'rbc')
            f.rbcshow = vin;
            varargout = {set(f,'rbc',vin)};
        elseif strcmp(idx,'ubc')
            f.ubcshow = vin;
            varargout = {set(f,'ubc',vin)};
        elseif strcmp(idx,'dbc')
            f.dbcshow = vin;
            varargout = {set(f,'dbc',vin)};
        elseif strcmp(idx,'op')
            varargout = {set(f,'op',vin)};
        elseif strcmp(idx,'domain')
            varargout = {set(f,'domain',vin)};
        elseif strcmp(idx,'coeffs')
            varargout = {set(f,'coeffs',vin)};
        elseif strcmp(idx,'xorder')
            varargout = {set(f,'xorder',vin)};   
        elseif strcmp(idx,'yorder')
            varargout = {set(f,'yorder',vin)};
        elseif strcmp(idx,'U')
            varargout = {set(f,'U',vin)};
        elseif strcmp(idx,'S')
            varargout = {set(f,'S',vin)};
        elseif strcmp(idx,'V')
            varargout = {set(f,'V',vin)};
        end
end
end