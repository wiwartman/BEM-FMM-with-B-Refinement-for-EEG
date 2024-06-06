for k = 1:length(Z)
    clf;
    z = Z(k);
    %  Display NIFTI slice
    I = round(z/d1d2d3(3) + N1N2N3(3)/2);
    S = VT(:, :, I)';      %   choose the Z cross-section
    S = S(:, end:-1:1);    
    image([-DimensionX/2 +DimensionX/2], [-DimensionY/2 +DimensionY/2], S, 'CDataMapping', 'scaled');
    colormap bone; brighten(0.3);
    set(gca, 'YDir', 'normal');    

    %   Create coordinates of intersection contours and intersection edges
    tissues = length(tissue);
    PofXY = cell(tissues, 1);   %   intersection nodes for a tissue
    EofXY = cell(tissues, 1);   %   edges formed by intersection nodes for a tissue
    TofXY = cell(tissues, 1);   %   intersected triangles
    NofXY = cell(tissues, 1);   %   normal vectors of intersected triangles
    count = [];   %   number of every tissue present in the slice
    for m = 1:tissues 
        [Pi, ti, polymask, flag] = meshplaneintXY(PS{m}, tS{m}, eS{m}, TriPS{m}, TriMS{m}, z);    
        if flag % intersection found                
            count               = [count m];
            PofXY{m}            = Pi;               %   intersection nodes
            EofXY{m}            = polymask;         %   edges formed by intersection nodes
            TofXY{m}            = ti;               %   intersected triangles
            NofXY{m}            = nS{m}(ti, :);     %   normal vectors of intersected triangles        
        end
    end
    %   Display the contours    
    for m = count
        edges           = EofXY{m};             %   this is for the contour
        points          = [];
        points(:, 1)    = +PofXY{m}(:, 1);       %   this is for the contour  
        points(:, 2)    = +PofXY{m}(:, 2);       %   this is for the contour
        ecenter         = (points(:, 1) + points(:, 2))/2;
        patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 3.0);    %   this is contour plot
        % hold on
        % plot(points(:, 1), points(:, 2), '.r', 'MarkerSize', 12);
    end
    %   Draw dipole
    hold on
    bemf1_graphics_dipole(1e3*strdipolePplus, 1e3*strdipolePminus, 1e3*strdipoleCurrent, 1);
    title( strcat('Transverse cross-section at z =', num2str(z), ' mm'));
    xlabel('x, mm'); ylabel('z, mm');
    axis 'equal';  axis 'tight'; 
    set(gcf,'Color','White');
    drawnow
    % pause(0.25)
 end