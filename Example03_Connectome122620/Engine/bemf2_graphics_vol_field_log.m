function bemf2_graphics_vol_field_log(temp, th1, th2, levels, a, b)
%   Volume field graphics:  plot a field quantity temp in the observation
%   plane using a planar contour plot with a log scale. Revision 05/03/24
%
%   temp - quantity to plot
%   th1, th2 - two threshold levels introduced manually
%   levels - number of levels in the contour plot introduced manually
%   a, b - x and y arguments
%
%   We apply the log-modulus transformation (John and Draper, 1980)
%   John JA, Draper NR. An Alternative Family of Transformations. 
%   J. of the Royal Statistical Society. Series C (Applied Statistics).
%   1980; 29(2): 190-197. doi: https://www.jstor.org/stable/2986305
%
%   Copyright SNM 2018-2024

    warning off %   to avoid warnings when the field is exactly zero
    temp(temp>+th1) = +th1;
    temp(temp<+th2) = +th2;     

    N = 11;
    factor          = 0.01; 
    scale           = factor*max(abs(temp));     
    templ           = sign(temp).*log10(abs(temp)/scale + 1);
    th1l            = sign(th1).*log10(abs(th1)/scale + 1);
    th2l            = sign(th2).*log10(abs(th2)/scale + 1);

    [C, h]          = contourf(a, b, reshape(templ, length(a), length(b)), levels);
    %tick            = round((th1l-th2l)/levels, 1, 'significant');
    % tick           = linspace(th1, th2, N);
    % h.LevelList     = scale*sign(tick).*(10.^(tick.*sign(tick))-1); %   this is the inverse functiontick*round(h.LevelList/tick);
    %h.LevelList     = tick*round(h.LevelList/tick);
    h.ShowText      = 'off';

    cb              = colorbar('FontSize', 15);
    cbscalein       = cb.Limits;
    cbscaleout      = [0 1];
    ticks           = linspace(cbscaleout(1),cbscaleout(2), N);
    cb.Ticks        = diff(cbscalein)*(ticks-cbscaleout(1))/diff(cbscaleout) + cbscalein(1);
    origvalues      = scale*sign(cb.Ticks).*(10.^(cb.Ticks.*sign(cb.Ticks))-1); %   this is the inverse function
    cb.TickLabels   = round(origvalues, 2, 'significant'); 

    warning on;

end