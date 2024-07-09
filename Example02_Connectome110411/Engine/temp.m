A = [10 20 30 20 50 30]

[~, uniqueIdx] = unique(A);
nondupeIdx = ~ismember(A, A(setdiff(1:numel(A), uniqueIdx)));