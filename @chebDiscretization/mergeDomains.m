function disc = mergeDomains(disc,varargin)

% NEWDSC = MERGEDOMAINS(DSC,D1,D2,...), where D1,... are numeric vectors,
% performs a chebfun merge of the numeric domains with that of the
% discretization, and assigns the result to the otherwise identical output
% domain NEWDSC.
%
% NEWDSC = MERGEDOMAINS(DSC,A1,A2,...) does the same using the domains of the
% chebmatrices A1,A2,....

if ( nargin == 1 )
    return
end

if isnumeric(varargin{1})
     d = varargin;
else
    d = { mergeDomains(varargin{:}) };
end

disc.domain = chebfun.mergeDomains(disc.domain,d{:});

end
