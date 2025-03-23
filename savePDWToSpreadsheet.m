function savePDWToSpreadsheet(data, baseFilename, fileExtension, folderPath)
    % Save a file with the current date and time appended to the filename in a specific folder
    % Inputs:
    %   data: The data to save (e.g., a table, matrix, etc.)
    %   baseFilename: The base name of the file (e.g., 'pdwTable')
    %   fileExtension: The file extension (e.g., '.xlsx', '.csv')
    %   folderPath: The path to the folder where the file should be saved

    % Get the current date and time
    currentTime = datetime('now', 'Format', 'yyyy_MMdd_HHmmss'); % Format: YearMonthDay_HourMinuteSecond
    
    % Convert the datetime to a string
    timestamp = char(currentTime);
    
    % Create the dynamic filename
    filename = sprintf('%s_%s%s', baseFilename, timestamp, fileExtension);
    
    % Construct the full file path
    fullFilePath = fullfile(folderPath, filename);
    
    % Save the data to the file
    if istable(data)
        writetable(data, fullFilePath); % For tables
    else
        writematrix(data, fullFilePath); % For matrices or arrays
    end
    
    % Display a confirmation message
    fprintf('File saved as %s\n', fullFilePath);
end