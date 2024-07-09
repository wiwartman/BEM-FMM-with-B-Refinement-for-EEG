function [I]= potint_self(r1, r2, r3)
%   Potential self integrals - surface patches -vectorized
%   r1 - Vertex1 (N, :), N - number of facets - vertex 1 of each triangle
%   r2 - Vertex1 (N, :), N - number of facets - vertex 2 of each triangle
%   r3 - Vertex1 (N, :), N - number of facets - vertex 3 of each triangle

%   I = IsIs(1/r)   1 integral total
%   T. F. Eibert and V. Hansen, 1995, "On the calculation of potential
%   integrals for linear source distributions on triangular domains," IEEE
%   Trans. Antennas Propagation, vol. 43, no. 12, Dec. 1995, pp. 1499-1502.
%   ECE539-ECE5106 April 2010 through May 2024 SNM ECE WPI, Worcester MA

r12 = r2-r1;
r23 = r3-r2;
r13 = r3-r1;

a = sum(r13.*r13, 2);       %   column, dot product
b = sum(r13.*r23, 2);       %   column, dot product
c = sum(r23.*r23, 2);       %   column, dot product
d = a-2*b+c;                %   column, dot product

A = sqrt(a);                %   column
B = sqrt(b);                %   column
C = sqrt(c);                %   column
D = sqrt(d);                %   column

%   Analytical formula for the first integral (Eibert, Hansen, 1995)
%   I(1/r)
N1 = (a-b+A.*D).*(b+A.*C);
D1 = (-a+b+A.*D).*(-b+A.*C);
    
N2 = (-b+c+C.*D).*(b+A.*C);
D2 = (b-c+C.*D).*(-b+A.*C);
    
N3 = (a-b+A.*D).*(-b+c+C.*D);
D3 = (b-c+C.*D).*(-a+b+A.*D);
    
Int     = 1/6*(1./A.*log(N1./D1) +1./C.*log(N2./D2) +1./D.*log(N3./D3));  
I       = 4 *Int; %  to obtain the full integral the result should be multiplied by A^2