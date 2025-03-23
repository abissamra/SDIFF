function [pulseTrain, pulseTrainMissing, pulseTrainWithNoise, pulseTrainMissingWithNoise] = generatePulseTrains(pulseConfig, fs, duration, noiseLevel)
    % generatePulseTrains - Gera trens de pulsos estáveis com e sem pulsos faltantes, com e sem ruído, e exibe os plots.
    %
    % Entradas:
    %   pulseConfig: Matriz de configurações dos pulsos. Cada linha representa [PRI, TP, Amplitude, Deslocamento, MissingPulses%].
    %   fs: Frequência de amostragem (Hz).
    %   duration: Duração total do sinal (s).
    %   noiseLevel: Nível do ruído (variância do ruído).
    %
    % Saídas:
    %   pulseTrain: Trem de pulsos entrelaçados sem pulsos faltantes e sem ruído.
    %   pulseTrainMissing: Trem de pulsos entrelaçados com pulsos faltantes e sem ruído.
    %   pulseTrainWithNoise: Trem de pulsos entrelaçados sem pulsos faltantes e com ruído.
    %   pulseTrainMissingWithNoise: Trem de pulsos entrelaçados com pulsos faltantes e com ruído.

    % Vetor de tempo
    t = 0:1/fs:duration;      

    % Inicializa os vetores de trem de pulsos
    pulseTrain = zeros(size(t));
    pulseTrainMissing = zeros(size(t));
    pulseTrainWithNoise = zeros(size(t));
    pulseTrainMissingWithNoise = zeros(size(t));
    numPulsesConfig = size(pulseConfig, 1);

    % Inicializa figuras
    figure; % Figura 1: Sem pulsos faltantes e sem ruído
    figure; % Figura 2: Com pulsos faltantes e sem ruído
    figure; % Figura 3: Sem pulsos faltantes e com ruído
    figure; % Figura 4: Com pulsos faltantes e com ruído

    % Gera e plota os trens de pulsos
    for j = 1:numPulsesConfig
        PRI = pulseConfig(j, 1);       % Intervalo de repetição do pulso
        TP = pulseConfig(j, 2);        % Duração do pulso
        A = pulseConfig(j, 3);         % Amplitude do pulso
        timeOffset = pulseConfig(j, 4); % Deslocamento específico para este pulso
        missingPulsesPercent = pulseConfig(j, 5); % Porcentagem de pulsos faltantes
        
        % Vetores do trem de pulsos individual
        individualPulse = zeros(size(t));
        individualPulseMissing = zeros(size(t));
        
        numPulses = floor((duration - timeOffset) / PRI); % Número de pulsos possíveis
        numMissingPulses = floor(numPulses * missingPulsesPercent / 100); % Número de pulsos faltantes

        % Gera uma lista de índices de pulsos faltantes (apenas uma vez)
        missingPulseIndices = randperm(numPulses, numMissingPulses);

        % Gera os pulsos no mesmo loop
        for i = 0:numPulses-1
            pulseStart = i * PRI + timeOffset;  % Início do pulso com deslocamento
            pulseEnd = pulseStart + TP;          % Fim do pulso
            pulseIndices = find(t >= pulseStart & t < pulseEnd);  % Índices do pulso

            % Adiciona o pulso ao trem de pulsos sem pulsos faltantes
            individualPulse(pulseIndices) = A;

            % Adiciona o pulso ao trem de pulsos com pulsos faltantes (se não estiver na lista de faltantes)
            if ~ismember(i, missingPulseIndices)
                individualPulseMissing(pulseIndices) = A;
            end
        end

        % Adiciona ruído aos pulsos individuais
        individualPulseWithNoise = individualPulse + noiseLevel * randn(size(t));
        individualPulseMissingWithNoise = individualPulseMissing + noiseLevel * randn(size(t));

        % Figura 1: Sem pulsos faltantes e sem ruído
        figure(1);
        subplot(numPulsesConfig + 1, 1, j);
        plot(t, individualPulse);
        title(['Pulso ' num2str(j) ' sem Pulsos Faltantes e sem Ruído']);
        xlabel('Tempo (s)');
        ylabel('Amplitude');
        grid on;
        axis([0 duration -0.5 1.5]);

        % Figura 2: Com pulsos faltantes e sem ruído
        figure(2);
        subplot(numPulsesConfig + 1, 1, j);
        plot(t, individualPulseMissing);
        title(['Pulso ' num2str(j) ' com Pulsos Faltantes e sem Ruído']);
        xlabel('Tempo (s)');
        ylabel('Amplitude');
        grid on;
        axis([0 duration -0.5 1.5]);

        % Figura 3: Sem pulsos faltantes e com ruído
        figure(3);
        subplot(numPulsesConfig + 1, 1, j);
        plot(t, individualPulseWithNoise);
        title(['Pulso ' num2str(j) ' sem Pulsos Faltantes e com Ruído']);
        xlabel('Tempo (s)');
        ylabel('Amplitude');
        grid on;
        axis([0 duration -0.5 1.5]);

        % Figura 4: Com pulsos faltantes e com ruído
        figure(4);
        subplot(numPulsesConfig + 1, 1, j);
        plot(t, individualPulseMissingWithNoise);
        title(['Pulso ' num2str(j) ' com Pulsos Faltantes e com Ruído']);
        xlabel('Tempo (s)');
        ylabel('Amplitude');
        grid on;
        axis([0 duration -0.5 1.5]);

        % Adiciona ao trem de pulsos geral sem pulsos faltantes
        pulseTrain = pulseTrain + individualPulse;  % Acumula para o último subplot
        pulseTrainWithNoise = pulseTrainWithNoise + individualPulseWithNoise;  % Acumula com ruído

        % Adiciona ao trem de pulsos geral com pulsos faltantes
        pulseTrainMissing = pulseTrainMissing + individualPulseMissing;  % Acumula para o último subplot
        pulseTrainMissingWithNoise = pulseTrainMissingWithNoise + individualPulseMissingWithNoise;  % Acumula com ruído
    end

    % Plota os trens de pulsos acumulados
    % Figura 1: Sem pulsos faltantes e sem ruído
    figure(1);
    subplot(numPulsesConfig + 1, 1, numPulsesConfig + 1);
    plot(t, pulseTrain);
    title('Trem de Pulsos Acumulado sem Pulsos Faltantes e sem Ruído');
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    grid on;
    axis([0 duration -0.5 1.5]);

    % Figura 2: Com pulsos faltantes e sem ruído
    figure(2);
    subplot(numPulsesConfig + 1, 1, numPulsesConfig + 1);
    plot(t, pulseTrainMissing);
    title('Trem de Pulsos Acumulado com Pulsos Faltantes e sem Ruído');
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    grid on;
    axis([0 duration -0.5 1.5]);

    % Figura 3: Sem pulsos faltantes e com ruído
    figure(3);
    subplot(numPulsesConfig + 1, 1, numPulsesConfig + 1);
    plot(t, pulseTrainWithNoise);
    title('Trem de Pulsos Acumulado sem Pulsos Faltantes e com Ruído');
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    grid on;
    axis([0 duration -0.5 1.5]);

    % Figura 4: Com pulsos faltantes e com ruído
    figure(4);
    subplot(numPulsesConfig + 1, 1, numPulsesConfig + 1);
    plot(t, pulseTrainMissingWithNoise);
    title('Trem de Pulsos Acumulado com Pulsos Faltantes e com Ruído');
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    grid on;
    axis([0 duration -0.5 1.5]);
end