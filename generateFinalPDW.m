function pdwTable = generateFinalPDW(pdwTable, priValues, toaValues)
    % Add a new column for PRI values to the PDW table
    pdwTable.PRI = NaN(height(pdwTable), 1); % Initialize PRI column with NaN

    % Iterate through each PRI sequence
    for i = 1:length(priValues)
        % Get the PRI value for the current sequence
        pri = priValues(i);
        
        % Get the TOA values for the current sequence
        toaSequence = toaValues{i};
        
        % Iterate through each TOA value in the sequence
        for j = 1:length(toaSequence)
            % Find the row in the PDW table where TOA matches
            matchIdx = find(abs(pdwTable.TOA - toaSequence(j)) < 1e-6); % Tolerance for floating-point comparison
            
            % If a match is found, add the PRI value to the corresponding row
            if ~isempty(matchIdx)
                pdwTable.PRI(matchIdx) = pri;
            end
        end
    end
end