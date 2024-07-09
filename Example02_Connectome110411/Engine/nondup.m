function [Id] = nondup(t)
%   Returns logical indexes for non-repeated elements in a row A
%   SNM 2024
    [~, uniqueIdx]  = unique(t);
    tt              = t(setdiff(1:length(t), uniqueIdx));
    Id = ~ismember(t, tt);
end