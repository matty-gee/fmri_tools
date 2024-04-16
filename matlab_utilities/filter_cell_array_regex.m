function match_indices = filter_cell_array_regex(cell_array, patterns)

    % filter a cell array based on matches to a set of paterns
    % from: https://www.mathworks.com/matlabcentral/answers/592129-how-to-filter-set-of-strings-with-set-of-patterns

    regex_matches = cell(numel(patterns), numel(cell_array));
    for i=1:numel(patterns)
        regex_matches(i, :) = regexpi(cell_array, patterns{i}, 'match');
    end
    match_indices = any(~cellfun(@isempty, regex_matches));