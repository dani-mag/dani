clc; close all;
load vep.mat;
Fs = 512;
Tventana = ((0:306) / Fs)*1000-300; 

%Para 120 ventanas
limite_ruido = 153; % 300 milisegundos (153 muestras)
Ventanas = zeros(120, 307, 10); % 120 ventanas, 307 muestras, 10 canales
promedios = zeros(10, 307); 

%Para 60 ventanas
Ventanas60=zeros(60, 307, 10);
promedios60 = zeros(10, 307);

%Para 20 ventanas
Ventanas20=zeros(20, 307, 10);
promedios20 = zeros(10, 307);

% 120 ventanas de los 10 canales
for i = 1:120
    for canal = 1:10
        Ventanas(i, :, canal) = obtenerVentana(canal, EEG, stim_times(i));
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

% Promedio de las 120 ventanas de los 10 canales
for canal = 1:10
    promedios(canal, :) = mean(Ventanas(:, :, canal), 1);
end

% Promedio de las primeras 60 ventanas de los 10 canales
for canal = 1:10
    promedios60(canal, :) = mean(Ventanas60(:, :, canal), 1);
end

% Promedio de las primeras 20 ventanas de los 10 canales
for canal = 1:10
    promedios20(canal, :) = mean(Ventanas20(:, :, canal), 1);
end

% Calcular el SNR para 120 ventanas
for canal = 1:10
    ruido = promedios(canal, 1:limite_ruido); % antes de 300 ms
    senal = promedios(canal, limite_ruido+1:end-1); % después de 300 ms
    SNR_value = snr(senal, ruido);

    % Graficar promedios con SNR
    figure();
    plot(Tventana, promedios(canal, :));
    title(['PEV para 120 ventanas del Canal ', num2str(canal)]);
    xlabel('Tiempo (s)');
    ylabel('Amplitud Promedio');
    text(Tventana(10), max(promedios(canal, :)) * 0.9, ['SNR: ', num2str(SNR_value, '%.2f'), ' dB'], 'FontSize', 12, 'BackgroundColor', 'w');

    % P100
    [pico, indice_pico] = max(promedios(canal, :)); 
    tiempo_pico = Tventana(indice_pico); 
    hold on;
    plot([tiempo_pico, tiempo_pico], ylim, '--', 'Color', [1, 0.2, 0.6], 'LineWidth', 1.5); 
    legend('PEV', 'P100', 'Location', 'northeast');
    grid on;
end

% Calcular el SNR para 60 ventanas
for canal = 1:10
    ruido60 = promedios60(canal, 1:limite_ruido); % antes de 300 ms
    senal60 = promedios60(canal, limite_ruido+1:end-1); % después de 300 ms
    SNR_value_60 = snr(senal60, ruido60);

    % Graficar promedios con SNR
    figure();
    plot(Tventana, promedios60(canal, :));
    title(['PEV para 60 ventanas del Canal ', num2str(canal)]);
    xlabel('Tiempo (s)');
    ylabel('Amplitud Promedio');
    text(Tventana(10), max(promedios60(canal, :)) * 0.9, ['SNR: ', num2str(SNR_value_60, '%.2f'), ' dB'], 'FontSize', 12, 'BackgroundColor', 'w');

    % P100
    [pico, indice_pico] = max(promedios60(canal, :)); 
    tiempo_pico = Tventana(indice_pico); 
    hold on;
    plot([tiempo_pico, tiempo_pico], ylim, '--', 'Color', [1, 0.2, 0.6], 'LineWidth', 1.5); 
    legend('PEV', 'P100', 'Location', 'northeast');
    grid on;
end

% Calcular el SNR para 20 ventanas
for canal = 1:10
    ruido20 = promedios20(canal, 1:limite_ruido); % antes de 300 ms
    senal20 = promedios20(canal, limite_ruido+1:end-1); % después de 300 ms
    SNR_value_20 = snr(senal20, ruido20);

    % Graficar promedios con SNR
    figure();
    plot(Tventana, promedios20(canal, :));
    title(['PEV para 20 ventanas del Canal ', num2str(canal)]);
    xlabel('Tiempo (s)');
    ylabel('Amplitud Promedio');
    text(Tventana(10), max(promedios20(canal, :)) * 0.9, ['SNR: ', num2str(SNR_value_20, '%.2f'), ' dB'], 'FontSize', 12, 'BackgroundColor', 'w');

    % P100
    [pico, indice_pico] = max(promedios20(canal, :)); 
    tiempo_pico = Tventana(indice_pico); 
    hold on;
    plot([tiempo_pico, tiempo_pico], ylim, '--', 'Color', [1, 0.2, 0.6], 'LineWidth', 1.5); 
    legend('PEV', 'P100', 'Location', 'northeast');
    grid on;
end

% Función para obtener las ventanas
function V = obtenerVentana(i, EEG, stim)
    V = EEG(i, (stim - 153):(stim + 153)); 
end