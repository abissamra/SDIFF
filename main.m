%Limpa variáveis e fecha abas abertas
close all
clear all
clc

% Parâmetros do trem de pulsos
fs = 1000;                  % Frequência de amostragem (Hz)
duration = 8;               % Duração total do sinal (s)
noiseLevel = 0.1;           % Nível do ruído (variância do ruído)

%Parâmetros de detecção dos pulsos
threshold = 0.1; % Limiar para detecção de pulsos
minPulseWidth = 0.00001; % Largura mínima do pulso (em segundos)

%Parâmetro de agrupamento por PW
pwTolerance = 0.002;

% Configurações dos pulsos
% Cada linha representa: [PRI, TP, Amplitude, Deslocamento, MissingPulses%]
pulseConfig = [
    0.11, 0.004, 1, 0.11, 10;  % Pulso 1: PRI = 0.6s, TP = 0.04s, Amplitude = 1, Deslocamento = 0.15s, MissingPulses% = 10%
    0.13, 0.004, 0.8, 0.13, 20; % Pulso 2: PRI = 0.8s, TP = 0.04s, Amplitude = 0.8, Deslocamento = 0.2s, MissingPulses% = 20%
    0.17, 0.003, 0.6, 0.17, 30;  % Pulso 3: PRI = 1.0s, TP = 0.03s, Amplitude = 0.6, Deslocamento = 0.3s, MissingPulses% = 30%
    0.19, 0.003, 0.5, 0.19, 10;
    %];
    0.23, 0.003, 0.5, 0.23, 10;
];

% Chama a função para gerar os trens de pulsos e exibir os plots
[pulseTrain, pulseTrainMissing, pulseTrainWithNoise, pulseTrainMissingWithNoise] = generatePulseTrains(pulseConfig, fs, duration, noiseLevel);

% Chama a função para detecção de TOA e PW dos pulsos
[TOA, PW, PDW] = detectPulsesCalculatePDW(pulseTrain, fs, threshold, minPulseWidth);

TOA_reorganized = reorganizeTOA(TOA, PW, pwTolerance);

TOA_copy = TOA;

[priValues, toaValues] = sdiff_2(TOA_copy, fs, duration);

PDW = generateFinalPDW(PDW, priValues, toaValues);

% Display the PDW
disp('PDW');
disp(PDW);

% Specify the folder path
folderPath = 'C:\Users\pedro\OneDrive\Documentos\GitHub\SDIF\Analises';

%savePDWToSpreadsheet(PDW,'Analise_RADAR_Aquisicao_', '.xlsx', folderPath);

