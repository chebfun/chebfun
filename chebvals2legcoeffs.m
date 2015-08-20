function c_leg = chebvals2legcoeffs(v_cheb, varargin)

c_cheb = chebvals2chebcoeffs(v_cheb);
c_leg = cheb2leg(c_cheb, varargin{:});

end