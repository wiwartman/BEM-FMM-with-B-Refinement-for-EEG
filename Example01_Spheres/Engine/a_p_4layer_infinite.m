function Pot = a_p_4layer_infinite(I0, Pplus, Pminus, radfactor, cond, Points)
%   Analytical solution  for a four-layer sphere with isostropic
%   conductivities. Solution follows [Z. Zhang, �A fast method to compute
%   surface potentials generated by dipoles within multilayer anisotropic
%   spheres,� Phys. Med. Biol., vol. 40, pp. 335�349, Mar. 1995]. Note: In
%   the review [Mosher JC, Leahy RM, Lewis PS. EEG and MEG: forward
%   solutions for inverse methods. IEEE Trans Biomed Eng. 1999
%   Mar;46(3):245-59] a similar result is given, but in Eq. (17) it should
%   be r instead of rq.
%
%   Copyright SNM 2018-2020

% Inputs:
% I0        - current amplitude
% PPlus     - position of +I0 source
% PMinus    - position of -I0 source
% radfactor - sphere radii
% cond      - sphere conductivities, S/m
% Points    - Mx3 array of positions on the sphere where the potential is plotted
%   Output(s):
% Pot - surface electric potential

    %   Sphere data
    r     = 0.1*radfactor(end:-1:1);    % sphere radii in ascending order  vs. 100 mm 
    sig   = cond(end:-1:1);             % conductivities (innermost first)    
    M     = length(r);
    sigM  = sig(end);

    N = 100; %  Series length
    K = size(Points, 1); 
    Pot = zeros(K, 1);
    dvector     = Pplus - Pminus;
    dcenter     = 0.5*(Pplus + Pminus);
    if norm(cross(dvector, dcenter))<1e-9*norm(dvector)*norm(dcenter)
        dvector = dvector + 1e-12*norm(dvector)*[1 1 1];
    end
    alpha       = acos(dot(dcenter, dvector)/norm(dcenter)/norm(dvector));
    cosalpha    = cos(alpha);
    sinalpha    = sin(alpha);
    qrqnormal   = cross(dcenter, dvector);
    qrqnormal   = qrqnormal/norm(qrqnormal);
    q           = I0*norm(dvector);
    rqs         = norm(dcenter);

    f = zeros(N, 1);
    for n = 1:N  
        MATRIX3 = [n+(n+1)*sig(3)/sig(4)                     (n+1)*(sig(3)/sig(4)-1)*(r(end)/r(3))^(2*n+1);...
                   n*(sig(3)/sig(4)-1)*(r(3)/r(end))^(2*n+1) (n+1)+n*sig(3)/sig(4)];

        MATRIX2 = [n+(n+1)*sig(2)/sig(3)                     (n+1)*(sig(2)/sig(3)-1)*(r(end)/r(2))^(2*n+1);...
                   n*(sig(2)/sig(3)-1)*(r(2)/r(end))^(2*n+1) (n+1)+n*sig(2)/sig(3)];

        MATRIX1 = [n+(n+1)*sig(1)/sig(2)                     (n+1)*(sig(1)/sig(2)-1)*(r(end)/r(1))^(2*n+1);...
                   n*(sig(1)/sig(2)-1)*(r(1)/r(end))^(2*n+1) (n+1)+n*sig(1)/sig(2)];

        MATRIX = (1/(2*n+1)^(M-1))*MATRIX1*MATRIX2*MATRIX3;
        m22    = MATRIX(2, 2);
        m21    = MATRIX(2, 1);
        f(n)     = n/(n*m22 + (1+n)*m21); 
    end    
    COSGAMMA = zeros(1, K);
    L1       = zeros(N, K);
    L2       = zeros(N, K);
    for m = 1:K
        point     = Points(m, :);  %   Single observation point
        rs        = norm(point);
        COSGAMMA(m)  = dot(dcenter, point)/(rqs*rs);
    end
    for n = 1:N
        Legendre = legendre(n, COSGAMMA);
        L1(n, :) = Legendre(1, :);
        L2(n, :) = Legendre(2, :);
    end  
    for m = 1:K
        point     = Points(m, :);  %   Single observation point
        rs        = norm(point);        
        qrrnormal = cross(dcenter, point);
        qrrnormal = qrrnormal/norm(qrrnormal);
        cosbeta   = cos(pi - acos(dot(qrqnormal, qrrnormal)));
        Pot(m) = 0;
        factor = q/(4*pi*sigM*rs^2);
        for n = 1:N           
            Pot(m) = Pot(m) + factor*((2*n+1)/n)*(rqs/rs)^(n-1)*...
                f(n)*(n*cosalpha*L1(n, m) + cosbeta*sinalpha*L2(n, m));
        end
    end    
end