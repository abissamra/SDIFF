function [priValues, toaValues] = sdiff(TOA, fs, duration)
    % Parameters
    t = 0:1/fs:duration; % Time vector
    N = length(t)-1; % Maximum interval to consider
    P = length(TOA); % Total number of pulses
    thresh5 = 0.0007; % Threshold for sequence detection
    jitter_tolerance = 0.002; % Tolerance for grouping jittered PRIs
    
    % Constants for the new threshold function
    x = 0.04; % Experimentally determined constant
    k = 1; % Experimentally determined constant
    
    % Initialize output variables
    priValues = []; % Vector to store PRI values
    toaValues = {}; % Cell array to store TOA vectors for each PRI
    
    diff_level = 1; % Start with difference level 1
    pulses_left = P; % Initialize remaining pulses

    % Main loop: continue while more than 5 pulses remain
    while pulses_left > 5
        % Reset sequence_found flag for the current diff_level
        sequence_found = false;
    
        % Compute histogram for current diff_level (SDIF)
        HIST = zeros(1, N); % Reset histogram for each diff_level
        for pulse = 1:(P - diff_level)
            interval = (TOA(pulse + diff_level) - TOA(pulse)); % Compute time difference
            interval_Hist = int32(interval * 1000); % Convert to integer index
            if interval_Hist > 0 && interval_Hist <= N
                HIST(interval_Hist) = HIST(interval_Hist) + 1; % Update histogram
            end
        end
    
        % Find significant intervals exceeding the threshold
        significant_intervals = [];
        for interval_Hist = 1:round(N/5)
            % Calculate the new threshold iteratively
            E = P; % Total number of pulses
            c = diff_level; % Current difference level
            t = interval_Hist; % Current bin number
            threshold_hist = x * (E - c) * exp(-t / (k * N)); % New threshold function
            
            if HIST(interval_Hist) > threshold_hist
                significant_intervals = [significant_intervals, interval_Hist / 1000]; % Store significant intervals
            end
        end
    
        % Group similar PRI values within jitter tolerance
        grouped_intervals = group_intervals(significant_intervals, jitter_tolerance);
    
        % If no groups are formed (grouped_intervals is empty), increment diff_level and continue
        if isempty(grouped_intervals)
            diff_level = diff_level + 1;
            pulses_left = P - diff_level;
            continue; % Skip to the next iteration of the main loop
        end
        
        central_intervals = cellfun(@median, grouped_intervals);
        
        % Perform subharmonic check
        potential_intervals = subharmonic_check(HIST, central_intervals, threshold_hist, N);
                
    
        % Apply the rule based on diff_level
        if diff_level == 1
            % For diff_level = 1, check if there is only one group
            if length(potential_intervals) == 1
                
                % Perform sequence search for the potential interval
                [TOA, sequence_found, pri, toa_sequence] = sequence_simple(TOA, potential_intervals, P, N, thresh5);
                if sequence_found
                    % Store the PRI and TOA sequence
                    priValues = [priValues; pri]; % Append PRI value
                    toaValues{end+1} = toa_sequence; % Append TOA sequence
                    
                    % Update P and pulses_left after removing the sequence
                    P = length(TOA);
                    diff_level = 1;
                    pulses_left = P - diff_level;
                end
            else
                % More than one group: increment diff_level
                diff_level = diff_level + 1;
                pulses_left = P - diff_level;
                continue;
            end
        else
            % For diff_level > 1, perform sequence search for all significant values
            for possiblePRI = potential_intervals
                
                % Perform sequence search for the potential interval
                [TOA, sequence_found, pri, toa_sequence] = sequence_simple(TOA, possiblePRI, P, N, thresh5);
                if sequence_found
                    % Store the PRI and TOA sequence
                    priValues = [priValues; pri]; % Append PRI value
                    toaValues{end+1} = toa_sequence; % Append TOA sequence
                    
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
    % If the input array is empty, return an empty cell array
    if isempty(intervals)
        grouped_intervals = {};
        return;
    end
    
    % Sort intervals
    intervals = sort(intervals);
    
    % Initialize groups
    grouped_intervals = {};
    current_group = [intervals(1)];
    
    % Group intervals within tolerance
    for i = 2:length(intervals)
        if abs(intervals(i) - current_group(end)) <= tolerance %* current_group(end)
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
function potential_interval = subharmonic_check(HIST, central_intervals, threshold_hist, N)
    % Parameters
    tolerance = 0.001;  % Time tolerance for harmonic checks (1 µs)
    
    % Step 1: Find the histogram maximum PRI (strongest peak)
    [max_count, max_idx] = max(HIST);
    max_interval = max_idx / 1000;  % Convert bin index to seconds
    
    % Step 2: Check if the strongest peak exceeds the threshold
    if max_count > threshold_hist * N / max_idx
        % Case 1: Strong peak → treat as the true PRI
        potential_interval = max_interval;
        is_subharmonic = false;
    else
        % Case 2: Weak peak → proceed with subharmonic checks
        if isempty(central_intervals)
            % No valid intervals → return empty (fallback)
            potential_interval = [];
            is_subharmonic = false;
        else
            % Step 3: Check the LOWEST central interval first
            [lowest_interval, idx] = min(central_intervals);
            
            % Is it a harmonic (integer multiple) of the histogram max?
            remainder = mod(lowest_interval, max_interval);
            is_harmonic = (remainder < tolerance) || (abs(remainder - max_interval) < tolerance);
            
            if is_harmonic
                % Case 3a: Lowest interval is a harmonic → use it
                potential_interval = lowest_interval;
                is_subharmonic = true;
            else
                % Case 3b: Not a harmonic → return ALL central intervals for testing
                potential_interval = central_intervals;
                is_subharmonic = false;
            end
        end
    end
end

% ==================================================
% SEQUENCE Procedure: Simple Weighting Scheme
% ==================================================
function [TOA_updated, sequence_found, pri, toa_sequence] = sequence_simple(TOA, interval, P, N, threshold_cnt)
    % Initialize variables
    pulse_A = 1;
    pulse_B = 2;
    CNT = 0;
    sequence_found = false;
    pri = 0; % Initialize PRI
    toa_sequence = []; % Initialize TOA sequence

    % Define a small tolerance value for floating-point comparisons
    tolerance = 3e-3;

    % Arrays to store the indices/values of pulses in the sequence
    sequence_pulses = [];

    % Outer loop to restart the search for pulse_A and pulse_B
    while pulse_B <= P
        % Find first two pulses matching the interval
        while pulse_B <= P
            TOA_diff = TOA(pulse_B) - TOA(pulse_A);

            if abs(TOA_diff - interval) < tolerance
                % Record pulses A and B
                sequence_pulses = [sequence_pulses, pulse_A, pulse_B];
                break; % Found matching interval
            elseif TOA_diff < interval
                pulse_B = pulse_B + 1; % Move pulse_B forward
            else
                pulse_A = pulse_A + 1; % Move pulse_A forward
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
            TOA_diff = TOA(pulse_C) - TOA(pulse_B);

            if abs(TOA_diff - interval) < tolerance
                % Record pulse C
                sequence_pulses = [sequence_pulses, pulse_C];
                break; % Found matching interval
            elseif TOA_diff < interval
                pulse_C = pulse_C + 1; % Move pulse_C forward
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
        TOA_diff = TOA(pulse_D) - last_TOA;

        if abs(TOA_diff - interval) < tolerance
            % Record pulse D
            sequence_pulses = [sequence_pulses, pulse_D];
            % Update sequence
            PRI = PRI + TOA_diff;
            CNT = CNT + 1;
            last_TOA = TOA(pulse_D);
        elseif TOA_diff < interval
            pulse_D = pulse_D + 1; % Move pulse_D forward
        else
            last_TOA = last_TOA + PRI / CNT; % Extrapolate based on average PRI
        end
    end

    % Count of pulse pairs separated by the searched interval
    N_pairs = 0;
    for i = 1:length(sequence_pulses)-1
        for j = i+1:length(sequence_pulses)
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
        
        % Store the PRI and TOA sequence
        pri = interval; % PRI value
        toa_sequence = TOA(sequence_pulses); % TOA values for the sequence
        
        % Remove the pulses in the sequence from the TOA array
        TOA_updated = TOA(setdiff(1:P, sequence_pulses));
        sequence_found = true;
    else
        % Discard the recorded pulses if no valid sequence is found
        TOA_updated = TOA;
    end
end