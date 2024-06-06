function [Epri, Ppri] = bemf3_inc_field_electric_plain_dipoles(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, Points, prec)
%   Computes potential and electric field from the dipole distribution via the FMM
%   at observation points (Points)
%

    %   Define source (pole) positions and FMM pseudo charges         
    Positions   = 0.5*(strdipolePplus+strdipolePminus);
    d           = (strdipolePplus-strdipolePminus);
    I0oversigma = strdipoleCurrent(1:2:end)./strdipolesig(1:2:end);
    PseudoM     = +repmat(I0oversigma, 1, 3).*d;  % plus here!

    %   FMM 2019
    srcinfo.nd      = 1;                    %   one vector of charges  
    srcinfo.sources = Positions';           %   source points
    targ            = Points';              %   target points
    pg      = 0;                            %   nothing is evaluated at sources
    pgt     = 2;                            %   potential/field are evaluated at targets
    srcinfo.dipoles          = PseudoM.';   %   pseudo dipoles    
    U                        = lfmm3d(prec, srcinfo, pg, targ, pgt);
    Ppri                     = +1/(4*pi)*U.pottarg.';
    Epri(:, 1)               = -1/(4*pi)*U.gradtarg(1, :);
    Epri(:, 2)               = -1/(4*pi)*U.gradtarg(2, :);
    Epri(:, 3)               = -1/(4*pi)*U.gradtarg(3, :); 
end