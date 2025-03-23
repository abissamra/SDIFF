% SDIF / SEQUENCE Algorithm for Radar Pulse Deinterleaving
% ==================================================
% SDIF Procedure: Sequential Difference Histogram with Subharmonic Check
% ==================================================
function sdiff(TOA, fs, duration)
    % Parameters
    t = 0:1/fs:duration; % Time vector (from 0 to duration, with steps of 1/fs)
    N = length(t)-1; % Maximum interval to consider (number of bins in the histogram)
    P = length(TOA); % Total number of pulses in the input TOA array
    thresh4 = 0.4; % Threshold for histogram detection (used to identify significant intervals)
    thresh5 = 0.0007; % Threshold for sequence detection (used to validate PRI sequences)
    jitter_tolerance = 0.01; % Tolerance for grouping jittered PRIs (1% of PRI)
    
    diff_level = 1; % Start with difference level 1 (difference between consecutive pulses)
    pulses_left = P; % Initialize remaining pulses (starts with all pulses)

    % Main loop: continue while more than 5 pulses remain
    while P > 5
        % Reset sequence_found flag for the current diff_level
        sequence_found = false;

        % Adjust the threshold based on the number of remaining pulses
        x = 1-(P/N); % Scaling factor for the threshold
        threshold_hist = thresh4*x; % Adjusted threshold for histogram detection
        threshold_cnt = thresh5*x; % Adjusted threshold for sequence detection
    
        % Compute histogram for current diff_level (SDIF)
        HIST = zeros(1, N); % Initialize histogram with zeros
        for pulse = 1:(P - diff_level)
            % Compute time difference between pulses separated by diff_level
            interval = (TOA(pulse + diff_level) - TOA(pulse));
            
            % Convert the interval to an integer index for the histogram
            interval_Hist = int32(interval * 1000); % Multiply by 1000 to convert to milliseconds
            
            % Update the histogram if the interval is within valid range
            if interval_Hist > 0 && interval_Hist <= N
                HIST(interval_Hist) = HIST(interval_Hist) + 1; % Increment the histogram bin
            end
        end
    
        % Find significant intervals exceeding the threshold
        significant_intervals = []; % Initialize array to store significant intervals
        for interval_Hist = 1:round(N/5)
            % Check if the histogram value exceeds the threshold
            if HIST(interval_Hist) > threshold_hist * N / interval_Hist
                % Store the significant interval (convert back to seconds)
                significant_intervals = [significant_intervals, interval_Hist / 1000];
            end
        end
    
        % Group similar PRI values within jitter tolerance
        grouped_intervals = group_intervals(significant_intervals, jitter_tolerance);
    
        % Apply the rule based on diff_level
        if diff_level == 1
            % For diff_level = 1, check if there is only one group
            if length(grouped_intervals) == 1
                % Use the central value of the group as the potential PRI
                central_interval = mean(grouped_intervals{1});
                
                % Perform subharmonic check
                [potential_interval, is_subharmonic] = subharmonic_check(HIST, central_interval, threshold_hist, N);
                
                % Perform sequence search for the potential interval
                [TOA, sequence_found] = sequence_simple(TOA, potential_interval, P, N, threshold_cnt);
                if sequence_found
                    % Update P and pulses_left after removing the sequence
                    P = length(TOA);
                    diff_level = 1;
                    pulses_left = P - diff_level;
                end
            else
                % More than one group: increment diff_level
                diff_level = diff_level + 1;
                pulses_left = P - diff_level;
            end
        else
            % For diff_level > 1, perform sequence search for all significant values
            for group = grouped_intervals
                % Use the central value of the group as the potential jittered PRI
                central_interval = mean(group);
                
                % Perform subharmonic check
                [potential_interval, is_subharmonic] = subharmonic_check(HIST, central_interval, threshold_hist, N);
                
                % Perform sequence search for the potential interval
                [TOA, sequence_found] = sequence_simple(TOA, potential_interval, P, N, threshold_cnt);
                if sequence_found
                    % Update P and pulses_left after removing the sequence
                    P = length(TOA);
                    diff_level = 1;
                    pulses_left = P - diff_level;
                    break; % Exit the loop to recompute HIST for the updated TOA
                end
            end
        end
    
        % If no sequences were found after checking all groups, increment diff_level
        if ~sequence_found
            diff_level = diff_level + 1;
            pulses_left = P - diff_level;
        end
    end
end

% ==================================================
% Group Similar Intervals Within Jitter Tolerance
% ==================================================
function grouped_intervals = group_intervals(intervals, tolerance)
    % Sort intervals
    intervals = sort(intervals);
    
    % Initialize groups
    grouped_intervals = {}; % Cell array to store groups of intervals
    current_group = [intervals(1)]; % Start the first group with the first interval
    
    % Group intervals within tolerance
    for i = 2:length(intervals)
        % Check if the current interval is within tolerance of the last interval in the group
        if abs(intervals(i) - current_group(end)) <= tolerance * current_group(end)
            % Add to current group if within tolerance
            current_group = [current_group, intervals(i)];
        else
            % Start a new group
            grouped_intervals{end+1} = current_group;
            current_group = [intervals(i)];
        end
    end
    
    % Add the last group
    grouped_intervals{end+1} = current_group;
end

% ==================================================
% Subharmonic Check
% ==================================================
function [potential_interval, is_subharmonic] = subharmonic_check(HIST, interval, threshold_hist, N)
    % Find the histogram maximum
    [max_count, max_idx] = max(HIST);
    max_interval = max_idx / 1000; % Convert to seconds
    
    % Check if the maximum exceeds the threshold
    if max_count > threshold_hist * N / max_idx
        % Use the maximum interval as the potential PRI
        potential_interval = max_interval;
        is_subharmonic = true;
    else
        % Check if the input interval is a subharmonic of the histogram maximum
        if mod(interval, max_interval) == 0
            % It is a subharmonic, so use the histogram maximum as the potential PRI
            potential_interval = max_interval;
            is_subharmonic = true;
        else
            % It is not a subharmonic, so use the input interval
            potential_interval = interval;
            is_subharmonic = false;
        end
    end
end

% ==================================================
% SEQUENCE Procedure: Simple Weighting Scheme
% ==================================================
function [TOA_updated, sequence_found] = sequence_simple(TOA, interval, P, N, threshold_cnt)
    % Initialize variables
    pulse_A = 1; % Index of the first pulse in the sequence
    pulse_B = 2; % Index of the second pulse in the sequence
    CNT = 0; % Counter for the number of pulses in the sequence
    sequence_found = false; % Flag to indicate if a valid sequence is found

    % Define a small tolerance value for floating-point comparisons
    tolerance = 3e-3;

    % Arrays to store the indices/values of pulses in the sequence
    sequence_pulses = [];

    % Outer loop to restart the search for pulse_A and pulse_B
    while pulse_B <= P
        % Find first two pulses matching the interval
        while pulse_B <= P
            % Compute the time difference between pulse_B and pulse_A
            TOA_diff = TOA(pulse_B) - TOA(pulse_A);

            % Check if the difference matches the interval within tolerance
            if abs(TOA_diff - interval) < tolerance
                % Record pulses A and B
                sequence_pulses = [sequence_pulses, pulse_A, pulse_B];
                break; % Found matching interval
            elseif TOA_diff < interval
                % Move pulse_B forward if the difference is too small
                pulse_B = pulse_B + 1;
            else
                % Move pulse_A forward if the difference is too large
                pulse_A = pulse_A + 1;
                pulse_B = pulse_B + 1;
            end
        end

        % If no matching interval is found, exit the function
        if isempty(sequence_pulses)
            TOA_updated = TOA;
            return;
        end

        % Find third pulse matching the interval
        pulse_C = pulse_B + 1;
        PRI = TOA_diff; % Initialize PRI with the first interval

        while pulse_C <= P
            % Compute the time difference between pulse_C and pulse_B
            TOA_diff = TOA(pulse_C) - TOA(pulse_B);

            % Check if the difference matches the interval within tolerance
            if abs(TOA_diff - interval) < tolerance
                % Record pulse C
                sequence_pulses = [sequence_pulses, pulse_C];
                break; % Found matching interval
            elseif TOA_diff < interval
                % Move pulse_C forward if the difference is too small
                pulse_C = pulse_C + 1;
            else
                % Reset sequence and restart the search for pulse_A and pulse_B
                pulse_A = pulse_A + 1;
                pulse_B = pulse_B + 1;
                sequence_pulses = []; % Clear the recorded pulses
                break; % Go back to the outer loop to restart the search
            end
        end

        % If a valid sequence is found, exit the outer loop
        if length(sequence_pulses) >= 3
            break;
        end

        if pulse_C > length(TOA)
            TOA_updated = TOA;
            return;
        end
    end

    % If no valid sequence is found, return the original TOA array
    if length(sequence_pulses) < 3
        TOA_updated = TOA;
        return;
    end

    % Update PRI and counters
    PRI = PRI + TOA_diff;
    CNT = 2;
    pulse_D = pulse_C + 1;
    last_TOA = TOA(pulse_C);

    % Extrapolate sequence
    while pulse_D <= P
        % Compute the time difference between pulse_D and the last pulse in the sequence
        TOA_diff = TOA(pulse_D) - last_TOA;

        % Check if the difference matches the interval within tolerance
        if abs(TOA_diff - interval) < tolerance
            % Record pulse D
            sequence_pulses = [sequence_pulses, pulse_D];
            % Update sequence
            PRI = PRI + TOA_diff;
            CNT = CNT + 1;
            last_TOA = TOA(pulse_D);
        elseif TOA_diff < interval
            % Move pulse_D forward if the difference is too small
            pulse_D = pulse_D + 1;
        else
            % Extrapolate based on average PRI if the difference is too large
            last_TOA = last_TOA + PRI / CNT;
        end
    end

    % Count of pulse pairs separated by the searched interval
    N_pairs = 0;
    for i = 1:length(sequence_pulses)-1
        for j = i+1:length(sequence_pulses)
            % Compute the time difference between pulses in the sequence
            TOA_diff = TOA(sequence_pulses(j)) - TOA(sequence_pulses(i));
            if abs(TOA_diff - interval) < tolerance
                N_pairs = N_pairs + 1;
            end
        end
    end

    % Compute the weighting function
    W = CNT + N_pairs;

    % Check if sequence is significant
    if W > threshold_cnt*N/interval % Threshold for sequence significance
        disp('Radar source found!'); % Declare radar source found
        disp(['Interval: ', num2str(interval)]);
        
        % Remove the pulses in the sequence from the TOA array
        TOA_updated = TOA(setdiff(1:P, sequence_pulses));
        sequence_found = true;
    else
        % Discard the recorded pulses if no valid sequence is found
        TOA_updated = TOA;
    end
end