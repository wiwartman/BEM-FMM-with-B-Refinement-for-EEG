for k = 1:length(X)
    clf;
    x = X(k);
    %   Display NIFTI slice
    I = round(-x/d1d2d3(1) + N1N2N3(1)/2); %    minus here!
    S = squeeze(VT(I, :, :))';      %   choose the X cross-section
    S = S(:, :);
    image([-DimensionY/2 +DimensionY/2], [-DimensionZ/2 +DimensionZ/2], S, 'CDataMapping', 'scaled');
    colormap bone; brighten(0.3);
    set(gca, 'YDir', 'normal');
    %   Create coordinates of intersection contours and intersection edges
    tissues = length(tissue);
    PofYZ = cell(tissues, 1);   %   intersection nodes for a tissue
    EofYZ = cell(tissues, 1);   %   edges formed by intersection nodes for a tissue
    TofYZ = cell(tissues, 1);   %   intersected triangles
    NofYZ = cell(tissues, 1);   %   normal vectors of intersected triangles
    count = [];   %   number of every tissue present in the slice
    for m = 1:tissues
        [Pi, ti, polymask, flag] = meshplaneintYZ(PS{m}, tS{m}, eS{m}, TriPS{m}, TriMS{m}, x);
        if flag % intersection found                
            count               = [count m];
            PofYZ{m}            = Pi;               %   intersection nodes
            EofYZ{m}            = polymask;         %   edges formed by intersection nodes
            TofYZ{m}            = ti;               %   intersected triangles
            NofYZ{m}            = nS{m}(ti, :);     %   normal vectors of intersected triangles        
        end
    end
    %   Display the contours    
    for m = count
        edges           = EofYZ{m};              %   this is for the contour
        points          = [];
        points(:, 1)    = +PofYZ{m}(:, 2);       %   this is for the contour  
        points(:, 2)    = +PofYZ{m}(:, 3);       %   this is for the contour
        patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 3.0);    %   this is contour plot
        % hold on
        % plot(points(:, 1), points(:, 2), '.r', 'MarkerSize', 12);
    end
    %   Draw dipole
    hold on
    bemf1_graphics_dipole(1e3*strdipolePplus, 1e3*strdipolePminus, strdipoleCurrent, 3);
    title( strcat('Sagittal cross-section at x =', num2str(x), ' mm'));
    xlabel('y, mm'); ylabel('z, mm');
    axis 'equal';  axis 'tight'; 
    set(gcf,'Color','White');
    drawnow
    %pause(0.25)
end
