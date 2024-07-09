function t = meshreorient(P, t, normals) 
%   This function reorients triangles
%   SNM 2024
    r1 = P(t(:, 1), :);
    r2 = P(t(:, 2), :);
    r3 = P(t(:, 3), :);
    tempv   = cross(r2-r1, r3-r1);      %   definition (*)
    tempd   = sum(tempv.*normals, 2);   %   dot product
    index   = tempd<0;
    t(index, 2:3) = t(index,  3:-1:2);  %   rearrange vertices to have exactly the outer normal by definition (*)
end