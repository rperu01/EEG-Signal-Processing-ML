%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BME772 Project -P04- ML to Detect Mental Arithmetic Tasks Using EEGs  %
% Ronah Peru - 500982554                                                %
% Jananan Velvelicham - XXXXXXXXX                                       %
% Mathew Szymanowski - XXXXXXXXX                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STUDENT A - Preprocessing Signals

%% FOR GRABBING THE EEG SIGNALS BEFORE

% Define the sampling frequency in Hz
fs = 500;

numberList = cell(1, 35);  % Preallocate space for 35 numbers

% Loop from 0 to 35
for i = 0:35
    % Format the number with two digits and leading zeros
    numberList{i+1} = sprintf('%02d', i);  % '%02d' ensures two digits with leading zeros
end

for i = 1:35
    % Defining the file name = point to file in with directory
    filename = ['C:\Users\ronah\OneDrive\Documents\MATLAB\BME772Project\eeg-during-mental-arithmetic-tasks-1.0.0\Subject', numberList{i},'_1.edf'];
    
    % Loading the data
    loadedData = edfread(filename);
    
    % Rename the 'Record Time' column to 'Time' for easier access
    loadedData.Properties.VariableNames{1} = 'Time';
    % Get the list of all lead names (excluding the 'Time' column)
    leadNames = loadedData.Properties.VariableNames(2:end);
    numLeads = numel(leadNames);
    
    % Total number of samples
    numSamples = size(loadedData, 1) * fs; 
    
    % Generate a continuous time vector
    timeVector = (0:numSamples-1) / fs;

    % Creating Variable for P3 signal
    EEGP3_B{i} = cell2mat(loadedData.(leadNames{12}));
    
    % Creating Variable for P4 signal
    EEGP4_B{i} = cell2mat(loadedData.(leadNames{13}));
end

%% FOR GRABBING THE EEG SIGNALS DURING

for i = 1:35
    % Defining the file name = point to file in with directory
    filename2 = ['C:\Users\ronah\OneDrive\Documents\MATLAB\BME772Project\eeg-during-mental-arithmetic-tasks-1.0.0\Subject', numberList{i},'_2.edf'];
    
    % Loading the data
    loadedData2 = edfread(filename2);
    
    % Rename the 'Record Time' column to 'Time' for easier access
    loadedData2.Properties.VariableNames{1} = 'Time';
    % Get the list of all lead names (excluding the 'Time' column)
    leadNames2 = loadedData2.Properties.VariableNames(2:end);
    numLeads2 = numel(leadNames);
    
    % Total number of samples
    numSamples = size(loadedData2, 1) * fs; 
    
    % Generate a continuous time vector
    timeVector = (0:numSamples-1) / fs;
    
    % Creating Variable for P3 signal
    EEGP3_D{i} = cell2mat(loadedData2.(leadNames{12}));
    
    % Creating Variable for P4 signal
    EEGP4_D{i} = cell2mat(loadedData2.(leadNames{13}));
end

%% Isolating Bands -EEGP3 Before
fs = 500; 

% Define frequency bands (in Hz)
bands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
frequencies = [0.5 4;    % Delta: 0.5 - 4 Hz
               4   7;    % Theta: 4 - 8 Hz
               8  12;    % Alpha: 8 - 13 Hz
               12 16;    % Beta: 13 - 30 Hz
               30 100];  % Gamma: 30 - 100 Hz (fixed upper bound)

% Assuming you have a list of filenames (or create them in sequence)
fileNames = cell(1, 35);  % Create a cell array for 35 files
for i = 1:35
    fileNames{i} = sprintf('file%02d.edf', i);  % Example filenames: file01.edf, file02.edf, ...
end

% Initialize the filtered signals cell array to hold signals for all bands
filteredSignals_P3B = cell(1, length(bands));  % Create a cell for each band

% Loop through each frequency band and filter the data
for m = 1:length(bands)

    % Create cell array to store filtered signals for each band (35 signals per band)
    filteredSignals_P3B{m} = cell(1, 35);  % Assuming 35 EEG signals for each band

    % Define the band for this iteration
    lowFreq = frequencies(m, 1);
    highFreq = frequencies(m, 2);
    
    % Design bandpass filter (3rd order Butterworth filter)
    [b, a] = butter(3, [lowFreq highFreq] / (fs / 2), 'bandpass');
    
    % Loop over each EEG signal
    for j = 1:35
        % Apply the filter to the EEG signal for channel j
        filteredSignals_P3B{m}{j} = filtfilt(b, a, EEGP3_B{j});  % Apply filter along the time dimension (columns)
    end
    
    % Calculate power for each band
    bandPowers_P3B = zeros(1, 35);  % Initialize an array to hold power values for each EEG signal
    
    % Loop through each EEG signal to calculate power in the current band
    for j = 1:35
        bandPowers_P3B(j) = mean(filteredSignals_P3B{m}{j}.^2);  % Power is mean squared amplitude of the signal
    end

    % Create a cell array to hold file numbers and corresponding band information
    fileNumbers = (1:35)';  % Assuming 35 EEG signals, with file numbers 1 to 35

    % Create a table for this band, adding a column for the file number
    bandTable = table(fileNumbers, repmat(bands(m), 35, 1), bandPowers_P3B', 'VariableNames', {'EEGP3B,Number', 'Band', 'Power'});
    
    % Display the table for the current frequency band
    disp(bandTable);

    % UNCOMMENT SECTION TO SEE ORIGINAL AND FILTERED SIGNALS FOR P3 BEFORE
    % Plot original and filtered signals for each EEG signal

    % for j = 1:3  % Just plot the first 3 EEG signals as an example
    %     figure;
    % 
    %     % Plot original signal
    %     subplot(2, 1, 1);
    %     plot(timeVector, EEGP3_B{j});
    %     title(['Original EEGP3 Before Signal - ' bands{m} ' - ' fileNames{j}]);
    %     xlabel('Time (s)');
    %     ylabel('Amplitude');
    %     axis tight;
    % 
    %     % Plot filtered signal
    %     subplot(2, 1, 2);
    %     plot(timeVector, filteredSignals_P3B{m}{j});
    %     title(['Filtered EEGP3 Before Signal - ' bands{m} ' - ' fileNames{j}]);
    %     xlabel('Time (s)');
    %     ylabel('Amplitude');
    %     axis tight;
    % end
end

%% Isolating Bands -EEGP4 Before

% Initialize the filtered signals cell array to hold signals for all bands
filteredSignals_P4B = cell(1, length(bands));  % Create a cell for each band

% Loop through each frequency band and filter the data
for m = 1:length(bands)

    % Create cell array to store filtered signals for each band (35 signals per band)
    filteredSignals_P4B{m} = cell(1, 35);  % Assuming 35 EEG signals for each band

    % Define the band for this iteration
    lowFreq = frequencies(m, 1);
    highFreq = frequencies(m, 2);
    
    % Design bandpass filter (3rd order Butterworth filter)
    [b, a] = butter(3, [lowFreq highFreq] / (fs / 2), 'bandpass');
    
    % Loop over each EEG signal
    for j = 1:35
        % Apply the filter to the EEG signal for channel j
        filteredSignals_P4B{m}{j} = filtfilt(b, a, EEGP4_B{j});  % Apply filter along the time dimension (columns)
    end
    
    % Calculate power for each band
    bandPowers_P4B = zeros(1, 35);  % Initialize an array to hold power values for each EEG signal
    
    % Loop through each EEG signal to calculate power in the current band
    for j = 1:35
        bandPowers_P4B(j) = mean(filteredSignals_P4B{m}{j}.^2);  % Power is mean squared amplitude of the signal
    end

    % Create a cell array to hold file numbers and corresponding band information
    fileNumbers = (1:35)';  % Assuming 35 EEG signals, with file numbers 1 to 35

    % Create a table for this band, adding a column for the file number
    bandTable = table(fileNumbers, repmat(bands(m), 35, 1), bandPowers_P4B', 'VariableNames', {'EEGP4B,Number', 'Band', 'Power'});
    
    % Display the table for the current frequency band
    disp(bandTable);

    % UNCOMMENT SECTION TO SEE ORIGINAL AND FILTERED SIGNALS FOR P4 BEFORE
    %Plot original and filtered signals for each EEG signal

    % for j = 1:3  % Just plot the first 3 EEG signals as an example
    %     figure;
    % 
    %     % Plot original signal
    %     subplot(2, 1, 1);
    %     plot(timeVector, EEGP4_B{j});
    %     title(['Original EEGP4 Before Signal - ' bands{m} ' - ' fileNames{j}]);
    %     xlabel('Time (s)');
    %     ylabel('Amplitude');
    %     axis tight;
    % 
    %     % Plot filtered signal
    %     subplot(2, 1, 2);
    %     plot(timeVector, filteredSignals_P4B{m}{j});
    %     title(['Filtered EEGP4 Before Signal - ' bands{m} ' - ' fileNames{j}]);
    %     xlabel('Time (s)');
    %     ylabel('Amplitude');
    %     axis tight;
    % end
end

%% Isolating Bands -EEGP3 During

% Initialize the filtered signals cell array to hold signals for all bands
filteredSignals_P3D = cell(1, length(bands));  % Create a cell for each band

% Loop through each frequency band and filter the data
for m = 1:length(bands)

    % Create cell array to store filtered signals for each band (35 signals per band)
    filteredSignals_P3D{m} = cell(1, 35);  % Assuming 35 EEG signals for each band

    % Define the band for this iteration
    lowFreq = frequencies(m, 1);
    highFreq = frequencies(m, 2);
    
    % Design bandpass filter (3rd order Butterworth filter)
    [b, a] = butter(3, [lowFreq highFreq] / (fs / 2), 'bandpass');
    
    % Loop over each EEG signal
    for j = 1:35
        % Apply the filter to the EEG signal for channel j
        filteredSignals_P3D{m}{j} = filtfilt(b, a, EEGP3_D{j});  % Apply filter along the time dimension (columns)
    end
    
    % Calculate power for each band
    bandPowers_P3D = zeros(1, 35);  % Initialize an array to hold power values for each EEG signal
    
    % Loop through each EEG signal to calculate power in the current band
    for j = 1:35
        bandPowers_P3D(j) = mean(filteredSignals_P3D{m}{j}.^2);  % Power is mean squared amplitude of the signal
    end

    % Create a cell array to hold file numbers and corresponding band information
    fileNumbers = (1:35)';  % Assuming 35 EEG signals, with file numbers 1 to 35

    % Create a table for this band, adding a column for the file number
    bandTable = table(fileNumbers, repmat(bands(m), 35, 1), bandPowers_P3D', 'VariableNames', {'EEGP3D,Number', 'Band', 'Power'});
    
    % Display the table for the current frequency band
    disp(bandTable);

    % UNCOMMENT SECTION TO SEE ORIGINAL AND FILTERED SIGNALS FOR P3 DURING 
    %Plot original and filtered signals for each EEG signal

    % for j = 1:3  % Just plot the first 3 EEG signals as an example
    %     figure;
    % 
    %     % Plot original signal
    %     subplot(2, 1, 1);
    %     plot(timeVector, EEGP3_D{j});
    %     title(['Original EEGP3 During Signal - ' bands{m} ' - ' fileNames{j}]);
    %     xlabel('Time (s)');
    %     ylabel('Amplitude');
    %     axis tight;
    % 
    %     % Plot filtered signal
    %     subplot(2, 1, 2);
    %     plot(timeVector, filteredSignals_P3D{m}{j});
    %     title(['Filtered EEGP3 During Signal - ' bands{m} ' - ' fileNames{j}]);
    %     xlabel('Time (s)');
    %     ylabel('Amplitude');
    %     axis tight;
    % end
end

%% Isolating Bands -EEGP4 During

% Initialize the filtered signals cell array to hold signals for all bands
filteredSignals_P4D = cell(1, length(bands));  % Create a cell for each band

% Loop through each frequency band and filter the data
for m = 1:length(bands)

    % Create cell array to store filtered signals for each band (35 signals per band)
    filteredSignals_P4D{m} = cell(1, 35);  % Assuming 35 EEG signals for each band

    % Define the band for this iteration
    lowFreq = frequencies(m, 1);
    highFreq = frequencies(m, 2);
    
    % Design bandpass filter (3rd order Butterworth filter)
    [b, a] = butter(3, [lowFreq highFreq] / (fs / 2), 'bandpass');
    
    % Loop over each EEG signal
    for j = 1:35
        % Apply the filter to the EEG signal for channel j
        filteredSignals_P4D{m}{j} = filtfilt(b, a, EEGP4_D{j});  % Apply filter along the time dimension (columns)
    end
    
    % Calculate power for each band
    bandPowers_P4D = zeros(1, 35);  % Initialize an array to hold power values for each EEG signal
    
    % Loop through each EEG signal to calculate power in the current band
    for j = 1:35
        bandPowers_P4D(j) = mean(filteredSignals_P4D{m}{j}.^2);  % Power is mean squared amplitude of the signal
    end

    % Create a cell array to hold file numbers and corresponding band information
    fileNumbers = (1:35)';  % Assuming 35 EEG signals, with file numbers 1 to 35

    % Create a table for this band, adding a column for the file number
    bandTable = table(fileNumbers, repmat(bands(m), 35, 1), bandPowers_P4D', 'VariableNames', {'EEGP4D,Number', 'Band', 'Power'});
    
    % Display the table for the current frequency band
    disp(bandTable);

    % UNCOMMENT SECTION TO SEE ORIGINAL AND FILTERED SIGNALS FOR P4 DURING   
    %Plot original and filtered signals for each EEG signal
    % 
    % for j = 1:3  % Just plot the first 3 EEG signals as an example
    %     figure;
    % 
    %     % Plot original signal
    %     subplot(2, 1, 1);
    %     plot(timeVector, EEGP4_D{j});
    %     title(['Original EEGP4 During Signal - ' bands{m} ' - ' fileNames{j}]);
    %     xlabel('Time (s)');
    %     ylabel('Amplitude');
    %     axis tight;
    % 
    %     % Plot filtered signal
    %     subplot(2, 1, 2);
    %     plot(timeVector, filteredSignals_P4D{m}{j});
    %     title(['Filtered EEGP4 During Signal - ' bands{m} ' - ' fileNames{j}]);
    %     xlabel('Time (s)');
    %     ylabel('Amplitude');
    %     axis tight;
    % end
end

%% STUDENT B - Feature Extraction

%% STUDENT C - Machine Learning
