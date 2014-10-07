function bcOut = bcReform(dom, bcInput, isIorF)
%BCREFORM    Convert a general BC from chebgui to IVP/FVP form
%
% BCOUT = BCREFORM(DOM, BCINPUT, ISIORF), where
%
%   DOM:     The domain that the problem is specified.
%   BCINPUT: The input from the BC field of CHEBGUI.
%   ISIORF:  Is equal to 1 if we're solving an initial value problem, and 2
%            if we're solving a final value problem.
%
% returns
%   BCOUT:   Boundary conditions that can be assigned to the LBC or RBC field of
%            a chebop.
%
% This method is used for converting input that usually would be assigned to
% the BC field of a chebop to form that can be assigned to the LBC or RBC field.
% It is required since in CHEBGUI, we specify conditions for ODEs on the format 
%   u(0) = 1,
% which can then be imposed via 
%   N.bc = @(x,u) u(0) - 1;
% On the other hand, for IVPs/FVPs, the format
%   u = 1
% is more useful, since that can then be imposed via
%   N.lbc = @(u) u;

% Find where whitespace and commas appear in the domain string, so that we can
% determine what the endpoints of the domain are:
spaceLoc = union(strfind(dom, ' '), strfind(dom, ','));

% Look at the BCINPUT string, and throw away parts of the string that we need to
% get rid of so that we convert it to an anonymous function that is suitable for
% N.lbc/rbc. Essentially, this is converting a string such as 
%   u(0) = 1
% to
%   u = 1.
% We need to throw different parts of the BCINPUT string away, depending on
% whether we're solving an IVP or FVP:
if ( isIorF ==  1 )
    % Initial value problem.
    firstSpace = min(spaceLoc);
    leftDomStr = dom(2:firstSpace-1);
    % String that we're going to be throwing away:
    throwString = ['(' ,leftDomStr, ')'];
else
    % Final value problem.
    lastSpace = max(spaceLoc);
    rightDomStr = dom(lastSpace+1:end-1);
    throwString = ['(' ,rightDomStr, ')'];
end

% Clear out the parts of the string required!
bcOut = strrep(bcInput, throwString, '');

end