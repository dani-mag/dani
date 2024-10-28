clc; close all;
load vep.mat;
Fs = 512;
Tventana = ((0:306) / Fs)*1000-300;
limite_ruido = 153; % 300 milisegundos= 153 muestras

%Matriz para 120 ventanas
Ventanas120 = zeros(120, 307, 10); % 120 ventanas, 307 muestras, 10 canales
Promedios120 = zeros(10, 307); 

%Matriz para 60 ventanas
Ventanas60=zeros(60, 307, 10); % 60 ventanas, 307 muestras, 10 canales
Promedios60 = zeros(10, 307);

%Matriz para 20 ventanas
Ventanas20=zeros(20, 307, 10); % 20 ventanas, 307 muestras, 10 canales
Promedios20 = zeros(10, 307);

%% VENTANAS

% 120 ventanas de los 10 canales
for i = 1:120
    for canal = 1:10
        Ventanas120(i, :, canal) = obtenerVentana(canal, EEG, stim_times(i));
    end
end

% 60 ventanas de los 10 canales
for i = 1:60
    for canal = 1:10
        Ventanas60(i, :, canal) = obtenerVentana(canal, EEG, stim_times(i));
    end
end

% 20 ventanas de los 10 canales
for i = 1:20
    for canal = 1:10
        Ventanas20(i, :, canal) = obtenerVentana(canal, EEG, stim_times(i));
    end
end

%% PROMEDIOS PARA GRAFICAR PEV

% Promedio de las 120 ventanas de los 10 canales
for canal = 1:10
    Promedios120(canal, :) = mean(Ventanas120(:, :, canal), 1);
end

% Promedio de las primeras 60 ventanas de los 10 canales
for canal = 1:10
    Promedios60(canal, :) = mean(Ventanas60(:, :, canal), 1);
end

% Promedio de las primeras 20 ventanas de los 10 canales
for canal = 1:10
    Promedios20(canal, :) = mean(Ventanas20(:, :, canal), 1);
end

%% SNR, P100 Y N75

% Calcular el SNR para 120 ventanas
for canal = 1:10
    ruido120 = Promedios120(canal, 1:limite_ruido); % antes de 300 ms
    senal120 = Promedios120(canal, limite_ruido+1:end-1); % después de 300 ms
    SNR_value_120 = snr(senal120, ruido120);

    % Graficar promedios con SNR
    figure();
    plot(Tventana, Promedios120(canal, :));
    title(['PEV para 120 ventanas del Canal ', num2str(canal)]);
    xlabel('Tiempo (ms)');
    ylabel('Amplitud Promedio');
    text(Tventana(10), max(Promedios120(canal, :)) * 0.9, ['SNR: ', num2str(SNR_value_120, '%.2f'), ' dB'], 'FontSize', 12, 'BackgroundColor', 'w');

    % P100
    [pico, indice_pico] = max(Promedios120(canal, :)); 
    tiempo_pico = Tventana(indice_pico); 
    hold on;
    plot([tiempo_pico, tiempo_pico], ylim, '--', 'Color', [1, 0.2, 0.6], 'LineWidth', 1.5); 
    legend('PEV', 'P100', 'Location', 'northeast');
    grid on;

    % N75 
    ventana_previa = Promedios120(canal, 1:indice_pico - 1);
    [pico_negativo_120, indice_pico_negativo_120] = findpeaks(-ventana_previa);
    if ~isempty(pico_negativo_120)
        indice_pico_negativo_120 = indice_pico_negativo_120(end);
        tiempo_pico_negativo_120 = Tventana(indice_pico_negativo_120);
        hold on;
        plot([tiempo_pico_negativo_120, tiempo_pico_negativo_120], ylim, '--', 'Color', [0, 1, 0], 'LineWidth', 1.5); 
        legend('PEV', 'P100', 'N75', 'Location', 'northeast');
    end
    grid on;
end

% Calcular el SNR para 60 ventanas
for canal = 1:10
    ruido60 = Promedios60(canal, 1:limite_ruido); % antes de 300 ms
    senal60 = Promedios60(canal, limite_ruido+1:end-1); % después de 300 ms
    SNR_value_60 = snr(senal60, ruido60);

    % Graficar promedios con SNR
    figure();
    plot(Tventana, Promedios60(canal, :));
    title(['PEV para 60 ventanas del Canal ', num2str(canal)]);
    xlabel('Tiempo (ms)');
    ylabel('Amplitud Promedio');
    text(Tventana(10), max(Promedios60(canal, :)) * 0.9, ['SNR: ', num2str(SNR_value_60, '%.2f'), ' dB'], 'FontSize', 12, 'BackgroundColor', 'w');

    % P100
    [pico, indice_pico] = max(Promedios60(canal, :)); 
    tiempo_pico = Tventana(indice_pico); 
    hold on;
    plot([tiempo_pico, tiempo_pico], ylim, '--', 'Color', [1, 0.2, 0.6], 'LineWidth', 1.5); 
    legend('PEV', 'P100', 'Location', 'northeast');
    grid on;

    % N75
    ventana_previa_60 = Promedios60(canal, 1:indice_pico - 1); 
    [pico_negativo_60, indice_pico_negativo_60] = findpeaks(-ventana_previa_60); 
    if ~isempty(pico_negativo_60)
        indice_pico_negativo_60 = indice_pico_negativo_60(end); 
        tiempo_pico_negativo_60 = Tventana(indice_pico_negativo_60);
        hold on;
        plot([tiempo_pico_negativo_60, tiempo_pico_negativo_60], ylim, '--', 'Color', [0, 1, 0], 'LineWidth', 1.5); 
        legend('PEV', 'P100', 'N75', 'Location', 'northeast');
    end
    grid on;
end

% Calcular el SNR para 20 ventanas
for canal = 1:10
    ruido20 = Promedios20(canal, 1:limite_ruido); % antes de 300 ms
    senal20 = Promedios20(canal, limite_ruido+1:end-1); % después de 300 ms
    SNR_value_20 = snr(senal20, ruido20);

    % Graficar promedios con SNR
    figure();
    plot(Tventana, Promedios20(canal, :));
    title(['PEV para 20 ventanas del Canal ', num2str(canal)]);
    xlabel('Tiempo (ms)');
    ylabel('Amplitud Promedio');
    text(Tventana(10), max(Promedios20(canal, :)) * 0.9, ['SNR: ', num2str(SNR_value_20, '%.2f'), ' dB'], 'FontSize', 12, 'BackgroundColor', 'w');

    % P100
    [pico, indice_pico] = max(Promedios20(canal, :)); 
    tiempo_pico = Tventana(indice_pico); 
    hold on;
    plot([tiempo_pico, tiempo_pico], ylim, '--', 'Color', [1, 0.2, 0.6], 'LineWidth', 1.5); 
    legend('PEV', 'P100', 'Location', 'northeast');
    grid on;

    % N75
    ventana_previa_20 = Promedios20(canal, 1:indice_pico - 1); 
    [pico_negativo_20, indice_pico_negativo_20] = findpeaks(-ventana_previa_20); 
    if ~isempty(pico_negativo_20)
        indice_pico_negativo_20 = indice_pico_negativo_20(end); 
        tiempo_pico_negativo_20 = Tventana(indice_pico_negativo_20);
        hold on;
        plot([tiempo_pico_negativo_20, tiempo_pico_negativo_20], ylim, '--', 'Color', [0, 1, 0], 'LineWidth', 1.5); 
        legend('PEV', 'P100', 'N75', 'Location', 'northeast');
    end
    grid on;
end

% Función para obtener las ventanas
function V = obtenerVentana(i, EEG, stim)
    V = EEG(i, (stim - 153):(stim + 153)); 
end