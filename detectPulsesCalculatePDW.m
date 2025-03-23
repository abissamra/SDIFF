function [TOA, PW, pdwTable] = detectPulsesCalculatePDW(pulseFrame, fs, threshold, minPulseWidth)
    % detectPulsesAndCalculatePDW - Detecta pulsos e calcula o PDW (Pulse Description Word) para cada pulso.
    % Saídas: TOA, PW e pdwTable (com Pulse Number, TOA, PW, Amplitude).
    %
    % Entradas:
    %   pulseFrame: Vetor do frame de pulsos (sinal no domínio do tempo).
    %   fs: Frequência de amostragem (Hz).
    %   threshold: Limiar para detectar o início de um pulso.
    %   minPulseWidth: Largura mínima de um pulso para ser considerado (em segundos).
    %
    % Saídas:
    %   TOA: Vetor com os tempos de chegada de cada pulso (em segundos).
    %   PW: Vetor com as larguras de cada pulso (em segundos).
    %   pdwTable: Tabela contendo o PDW de cada pulso (Pulse Number, TOA, PW, Amplitude).

    % Vetor de tempo
    t = (0:length(pulseFrame)-1) / fs;

    % Inicializa variáveis
    TOA = [];
    PW = [];
    Amplitude = [];
    isPulseActive = false;
    pulseStartTime = 0;
    pulseNumber = 0;

    % Hysteresis thresholds
    thresholdStart = threshold; % Limiar para início do pulso
    thresholdEnd = threshold * 0.8; % Limiar para fim do pulso (menor que o limiar de início)

    % Percorre o frame para detectar os pulsos
    for i = 1:length(pulseFrame)
        if pulseFrame(i) > thresholdStart && ~isPulseActive
            % Início de um novo pulso
            isPulseActive = true;
            pulseStartTime = t(i);
            pulseAmplitude = pulseFrame(i); % Inicializa a amplitude do pulso
        elseif pulseFrame(i) > thresholdStart && isPulseActive
            % Atualiza a amplitude máxima do pulso
            pulseAmplitude = max(pulseAmplitude, pulseFrame(i));
        elseif pulseFrame(i) <= thresholdEnd && isPulseActive
            % Fim do pulso atual
            isPulseActive = false;
            pulseEndTime = t(i);
            pulseWidth = pulseEndTime - pulseStartTime;

            % Verifica se a largura do pulso é maior que o mínimo especificado
            if pulseWidth >= minPulseWidth
                pulseNumber = pulseNumber + 1; % Incrementa o número do pulso
                TOA = [TOA; pulseStartTime]; % Armazena o TOA
                PW = [PW; pulseWidth];       % Armazena a PW
                Amplitude = [Amplitude; pulseAmplitude]; % Armazena a amplitude
            end
        end
    end

    % Caso o último pulso não tenha terminado ao final do frame
    if isPulseActive
        pulseEndTime = t(end);
        pulseWidth = pulseEndTime - pulseStartTime;
        if pulseWidth >= minPulseWidth
            pulseNumber = pulseNumber + 1;
            TOA = [TOA; pulseStartTime];
            PW = [PW; pulseWidth];
            Amplitude = [Amplitude; pulseAmplitude];
        end
    end

    % Cria uma tabela para armazenar o PDW de cada pulso
    PulseNumber = (1:pulseNumber)';
    pdwTable = table(PulseNumber, TOA, PW, Amplitude, ...
        'VariableNames', {'PulseNumber', 'TOA', 'PW', 'Amplitude'});
end