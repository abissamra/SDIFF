function reorganized_toa = reorganizeTOA(TOA, PW, tolerance)
    % reorganizeTOA: Reorganizes TOA array into a matrix where each row contains
    % TOA values with similar PW.
    %
    % Inputs:
    %   TOA       - Array of Time of Arrival values.
    %   PW        - Array of Pulse Width values corresponding to TOA.
    %   tolerance - Tolerance for PW variation to group pulses.
    %
    % Output:
    %   reorganized_toa - Matrix where each row contains TOA values for pulses
    %                     with similar PW. Rows are padded with NaN if necessary.

    % Validate inputs
    if nargin < 3
        error('Not enough input arguments. TOA, PW, and tolerance are required.');
    end
    if length(TOA) ~= length(PW)
        error('TOA and PW must have the same length.');
    end

    % Group indices of pulses with similar PW
    grouped_indices = {};
    used_indices = false(1, length(PW));  % Track used indices

    for i = 1:length(PW)
        if ~used_indices(i)
            % Find all indices where PW is within the tolerance range
            similar_pw_indices = find(abs(PW - PW(i)) <= tolerance);
            grouped_indices{end+1} = similar_pw_indices;
            used_indices(similar_pw_indices) = true;  % Mark indices as used
        end
    end

    % Determine the maximum number of TOA values in any group
    max_group_size = max(cellfun(@length, grouped_indices));

    % Initialize the reorganized TOA matrix with NaN padding
    reorganized_toa = NaN(length(grouped_indices), max_group_size);

    % Fill the matrix with TOA values
    for i = 1:length(grouped_indices)
        group = TOA(grouped_indices{i});
        reorganized_toa(i, 1:length(group)) = group;
    end
end