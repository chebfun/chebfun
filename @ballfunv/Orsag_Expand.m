%% Apply the 3/2 rule
function Vexp = Orsag_Expand(V)
[m,n,p] = size(V);
m_tilde = ceil(3*m/2);
n_tilde = ceil(3*n/2);
p_tilde = ceil(3*p/2);

[Vx,Vy,Vz] = V.comp{:};
Vx = Vx.coeffs; Vy = Vy.coeffs; Vz = Vz.coeffs;
VexpX = zeros(m_tilde,n_tilde,p_tilde);
VexpX(1:m,floor(n_tilde/2)+1-floor(n/2):floor(n_tilde/2)+n-floor(n/2),floor(p_tilde/2)+1-floor(p/2):floor(p_tilde/2)+p-floor(p/2)) = Vx;
VexpY = zeros(m_tilde,n_tilde,p_tilde);
VexpY(1:m,floor(n_tilde/2)+1-floor(n/2):floor(n_tilde/2)+n-floor(n/2),floor(p_tilde/2)+1-floor(p/2):floor(p_tilde/2)+p-floor(p/2)) = Vy;
VexpZ = zeros(m_tilde,n_tilde,p_tilde);
VexpZ(1:m,floor(n_tilde/2)+1-floor(n/2):floor(n_tilde/2)+n-floor(n/2),floor(p_tilde/2)+1-floor(p/2):floor(p_tilde/2)+p-floor(p/2)) = Vz;
Vexp = ballfunv(ballfun(VexpX),ballfun(VexpY),ballfun(VexpZ));
end