classdef stringParser
%STRINGPARSER   Parse strings for CHEBGUI.
%   This class is not intended to be called directly by the end user.
%
%   See also CHEBGUI.
    
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    methods (Access = public)
        
        function A = stringParser(varargin)
            % We never construct objects of this type, so the constructor is
            % trivial
        end
        
    end
    
    methods( Access = public, Static = true )
        
        % The STRINGPARSER lexer
        [out, varNames, pdeVarNames, eigVarNames, indVarNames] = ...
            lexer(str, type)
    
        % The LL(1) parser
        parseOut = parser(lexIn)
        
        % Get rid of unecessary parenthesis in infix format strings
        strOut = parSimp(strIn)
        
        % Convert an expression on PREFIX form to INFIX form
        [infixOut, notaVAR] = pref2inf(prefIn)
    end
end