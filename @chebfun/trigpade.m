function varargout = trigpade(f, m, n, varargin)
%TRIGPADE   Trigonometric (Fourier) Pade approximation.
%   [P, Q, R_HANDLE] = TRIGPADE(F, M, N) computes trigonometric polynomials 
%   P and Q of degree M and N, respectively, such that the trigonometric 
%   rational function P/Q is the type (M,N) Fourier-Pade approximation of 
%   the periodic CHEBFUN F. That is, the Fourier series of P/Q coincides 
%   with that for the CHEBFUN F up to the maximum possible order for the 
%   polynomial degrees permitted. R_HANDLE is a function handle for 
%   evaluating the trigonometric rational function P/Q.
% 
%   [P, Q, R_HANDLE, TN_P, TD_P, TN_M, TD_M] = TRIGPADE(F, M, N) also
%   returns the four trigonometric polynomials TN_P, TD_P, TN_M and TD_M
%   such that P/Q = TN_P./TD_P + TN_M./TD_M.
%
%   In both of the above cases, if only one output argument is specified
%   then R_HANDLE is returned, while P and Q are returned if two or more 
%   output arguments are specified. 
%
%   References:
%
%   [1] G. A. Baker and P. R. Graves-Morris. Pade Approximants. 
%   Cambridge University Press, second edition, 1996.
%
%   [2] M. Javed. Algorithms for trigonometric polynomial and rational 
%   approximation. DPhil thesis, Oxford.
%
% See also PADEAPPROX and CHEBPADE

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if ( isempty(f) )
    varargout{1} = f;
    return
end

if ( ~isPeriodicTech(f) )
    error('CHEBFUN:CHEBFUN:trigpade:trig', ...
        'Input chebfun F must have a periodic representation');
end

dom = f.domain([1, end]);

% Extract the Fourier coefficients of F:
coeffs  = trigcoeffs(f);

% Index of the middle coefficint:
c = (length(coeffs)-1)/2;

% Separate the +ve and -ve coefficients
% _p is for 'plus' and _m is for minus
coeffs_p = coeffs(c+1:end);
coeffs_m = coeffs(c+1:-1:1);

% Half the zeroth order coefficient
coeffs_p(1) = 1/2*coeffs_p(1);
coeffs_m(1) = 1/2*coeffs_m(1);

% Solve the two Pade approximation problems
[~, a_p, b_p] = padeapprox(coeffs_p, m, n);
if ( isreal(f) )
    % If f is real then the second Pade approximant 
    % is the conjugate of the first.
    a_m = conj(a_p);
    b_m = conj(b_p);    
else 
    % If f is not real, then solve the second Pade problem
    % explicitly
    [~, a_m, b_m] = padeapprox(coeffs_m, m, n);           
end

%% Construct the four trigonometric polynomials of Fourier Pade approximaiton

% Pad coefficients with zeros:
aa_p = [zeros(length(a_p)-1, 1); a_p];
bb_p = [zeros(length(b_p)-1, 1); b_p];

% _d is for denominator, _n is for numerator
% Denonminator and numerator for the +ve part
tn_p = chebfun(aa_p, 'coeffs', 'trig' );
td_p = chebfun(bb_p, 'coeffs', 'trig' );

% Pad coefficients with zeros
aa_m = [zeros(length(a_m)-1, 1); a_m];
aa_m = flipud(aa_m);
bb_m = [zeros(length(b_m)-1, 1); b_m];
bb_m = flipud(bb_m);

% denonminator and numerator for the -ve part
tn_m = chebfun(aa_m, dom, 'coeffs', 'trig' );
td_m = chebfun(bb_m, dom, 'coeffs', 'trig' );

% Construct the full approximation:
p = tn_p.*td_m + tn_m.*td_p; 
q = td_p.*td_m;
r_h = @(t) p(t)./q(t);

% Discard the imaginary rounding errors:
if ( isreal(f) && (norm(imag(p)) + norm(imag(q)) > 100*norm(f)* eps ) )
	warning('CHEBFUN:CHEBFUN:trigpade:imag', 'imaginary part not negligible.');
else
	p = real(p);
    q = real(q);
    r_h = @(t) p(t)./q(t);
end

outArgs = {p, q, r_h, tn_p, td_p, tn_m, td_m};
if ( nargout <= 1 )
    varargout{1} = r_h;
elseif ( nargout <= 7 )
    [varargout{1:nargout}] = outArgs{1:nargout};
else
    error('CHEBFUN:CHEBFUN:trigpade:nargout', ...
        'Incorrect number of output arguments.'); 
end
    

end