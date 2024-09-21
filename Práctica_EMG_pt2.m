%% Cargar las señales y parámetros iniciales

signal_3lb = load('alex_3lb.mat');
signal_5lb = load('alex_5lb.mat');
signal_13lb = load('alex_13lb.mat');
signal_sostenido = load('alexsostenido_13lb.mat');

% Asignar las variables
emg_3lb = signal_3lb.data;  
emg_5lb = signal_5lb.data; 
emg_13lb = signal_13lb.data;  
emg_sostenido = signal_sostenido.data;  

% Extraer señales de bíceps (primera columna) y tríceps (segunda columna)
bicep_3lb = emg_3lb(:, 1);
tricep_3lb = emg_3lb(:, 2);

bicep_5lb = emg_5lb(:, 1);
tricep_5lb = emg_5lb(:, 2);

bicep_13lb = emg_13lb(:, 1);
tricep_13lb = emg_13lb(:, 2);

bicep_sostenido = emg_sostenido(:, 1);
tricep_sostenido = emg_sostenido(:, 2);

% Frecuencia de muestreo
fs = 5000;
dt = 1 / fs; 
t_3lb = (0:length(emg_3lb)-1) * dt; 
t_5lb = (0:length(emg_5lb)-1) * dt;
t_13lb = (0:length(emg_13lb)-1) * dt; 
t_sostenido = (0:length(emg_sostenido)-1) * dt; 

% Graficar las señales de EMG
figure;
    %biceps
subplot(2, 4, 1);
plot(t_3lb, bicep_3lb);
title('EMG Bicep 3 lb');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2, 4, 2);
plot(t_5lb, bicep_5lb);
title('EMG Bicep 5 lb');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2, 4, 3);
plot(t_13lb, bicep_13lb);
title('EMG Bicep 13 lb');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2, 4, 4);
plot(t_sostenido, bicep_sostenido);
title('EMG Bicep Sostenido 13 lb');
xlabel('Tiempo (s)');
ylabel('Amplitud');

    % tríceps
subplot(2, 4, 5);
plot(t_3lb, tricep_3lb);
title('EMG Tríceps 3 lb');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2, 4, 6);
plot(t_5lb, tricep_5lb);
title('EMG Tríceps 5 lb');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2, 4, 7);
plot(t_13lb, tricep_13lb);
title('EMG Tríceps 13 lb');
xlabel('Tiempo (s)');
ylabel('Amplitud');

subplot(2, 4, 8);
plot(t_sostenido, tricep_sostenido);
title('EMG Tríceps Sostenido 13 lb');
xlabel('Tiempo (s)');
ylabel('Amplitud');

sgtitle('Señales EMG de Bíceps y Tríceps');

%Segmentos definidos manualmente
segmentos_bicep_3lb = [3553 8116; 15502 18983; 26227 29568; 37094 39494; 48244 50502];
segmentos_bicep_5lb = [4941 8287; 16997 19416; 27158 29255; 38247 40586; 50102 52078];
segmentos_bicep_13lb= [2141 8822; 16395 20865; 27780 33049; 40576 45609; 53277 58264];

segmentos_tricep_3lb= [10374 13385; 21664 25004; 32672 34554; 42975 45609; 53983 56147];
segmentos_tricep_5lb= [11634 14578; 22602 24457; 32682 35021; 43811 47158; 55908 56802];
segmentos_tricep_13lb= [12020 14749; 24957 27215; 37189 40293; 49608 52430; 62309 64253];

% Segmentos para 3 lb
figure;
for i = 1:5
    subplot(2, 5, i);  
    plot(t_3lb(segmentos_bicep_3lb(i,1):segmentos_bicep_3lb(i,2)), ...
         bicep_3lb(segmentos_bicep_3lb(i,1):segmentos_bicep_3lb(i,2)));
    title(['Segmento ' num2str(i) ' Bíceps']);
    xlabel('Tiempo [s]');
    ylabel('Amplitud [V]');
    ylim([-0.25 0.25]);
    
    subplot(2, 5, i+5);  
    plot(t_3lb(segmentos_tricep_3lb(i,1):segmentos_tricep_3lb(i,2)), ...
         tricep_3lb(segmentos_tricep_3lb(i,1):segmentos_tricep_3lb(i,2)));
    title(['Segmento ' num2str(i) ' Tríceps']);
    xlabel('Tiempo [s]');
    ylabel('Amplitud [V]');
    ylim([-0.25 0.25]);

    sgtitle('Segmentos de Bíceps y Tríceps para 3lb');
end

%Segmentos para 5lb 
figure;
for i = 1:5
    subplot(2, 5, i);  
    plot(t_5lb(segmentos_bicep_5lb(i,1):segmentos_bicep_5lb(i,2)), ...
         bicep_5lb(segmentos_bicep_5lb(i,1):segmentos_bicep_5lb(i,2)));
    title(['Segmento ' num2str(i) ' Bíceps']);
    xlabel('Tiempo [s]');
    ylabel('Amplitud [V]');
    ylim([-0.25 0.25]);
    
    subplot(2, 5, i+5);  
    plot(t_5lb(segmentos_tricep_5lb(i,1):segmentos_tricep_5lb(i,2)), ...
         tricep_5lb(segmentos_tricep_5lb(i,1):segmentos_tricep_5lb(i,2)));
    title(['Segmento ' num2str(i) ' Tríceps']);
    xlabel('Tiempo [s]');
    ylabel('Amplitud [V]');
    ylim([-0.25 0.25]);

    sgtitle('Segmentos de Bíceps y Tríceps para 5lb');
end

%Segmentos para 13 lb
figure;
for i = 1:5
    subplot(2, 5, i);  
    plot(t_13lb(segmentos_bicep_13lb(i,1):segmentos_bicep_13lb(i,2)), ...
         bicep_13lb(segmentos_bicep_13lb(i,1):segmentos_bicep_13lb(i,2)));
    title(['Segmento ' num2str(i) ' Bíceps']);
    xlabel('Tiempo [s]');
    ylabel('Amplitud [V]');
    ylim([-0.25 0.25]);
    
    subplot(2, 5, i+5);  
    plot(t_13lb(segmentos_tricep_13lb(i,1):segmentos_tricep_13lb(i,2)), ...
         tricep_13lb(segmentos_tricep_13lb(i,1):segmentos_tricep_13lb(i,2)));
    title(['Segmento ' num2str(i) ' Tríceps']);
    xlabel('Tiempo [s]');
    ylabel('Amplitud [V]');
    ylim([-0.25 0.25]);

    sgtitle('Segmentos de Bíceps y Tríceps para 13lb');
end

%% CALCULOS DE RMS Y STD 

% Función para calcular el RMS
rms_valor = @(x) sqrt(mean(x.^2));

% Inicializar variables 
rms_bicep_3lb = zeros(1, 5);
rms_tricep_3lb = zeros(1, 5);
rms_bicep_5lb = zeros(1, 5);
rms_tricep_5lb = zeros(1, 5);
rms_bicep_13lb = zeros(1, 5);
rms_tricep_13lb = zeros(1, 5);

area_biceps_3lb = zeros(1, 5);
area_triceps_3lb = zeros(1, 5);
area_biceps_5lb = zeros(1, 5);
area_triceps_5lb = zeros(1, 5);
area_biceps_13lb = zeros(1, 5);
area_triceps_13lb = zeros(1, 5);

cruces_biceps_sostenido = zeros(1, 5);
cruces_triceps_sostenido = zeros(1, 5);

% Segmentos para 3 lb
disp('RMS para 3 lb:');
for i = 1:5
    % Bíceps
    segmento_bicep = bicep_3lb(segmentos_bicep_3lb(i,1):segmentos_bicep_3lb(i,2));
    rms_bicep_3lb(i) = rms_valor(segmento_bicep);
%     fprintf('RMS Bíceps Segmento %d: %f\n', i, rms_bicep_3lb(i));
    
    % Tríceps
    segmento_tricep = tricep_3lb(segmentos_tricep_3lb(i,1):segmentos_tricep_3lb(i,2));
    rms_tricep_3lb(i) = rms_valor(segmento_tricep);
%     fprintf('RMS Tríceps Segmento %d: %f\n', i, rms_tricep_3lb(i));
end

% Promedio y desviación estándar para 3 lb
mean_rms_bicep_3lb = mean(rms_bicep_3lb);
std_rms_bicep_3lb = std(rms_bicep_3lb);
mean_rms_tricep_3lb = mean(rms_tricep_3lb);
std_rms_tricep_3lb = std(rms_tricep_3lb);

fprintf('Promedio RMS Bíceps 3lb: %f, Desviación Estándar: %f\n', mean_rms_bicep_3lb, std_rms_bicep_3lb);
fprintf('Promedio RMS Tríceps 3lb: %f, Desviación Estándar: %f\n', mean_rms_tricep_3lb, std_rms_tricep_3lb);

% Segmentos para 5 lb
disp('RMS para 5 lb:');
for i = 1:5
    % Bíceps
    segmento_bicep = bicep_5lb(segmentos_bicep_5lb(i,1):segmentos_bicep_5lb(i,2));
    rms_bicep_5lb(i) = rms_valor(segmento_bicep);
%     fprintf('RMS Bíceps Segmento %d: %f\n', i, rms_bicep_5lb(i));
    
    % Tríceps
    segmento_tricep = tricep_5lb(segmentos_tricep_5lb(i,1):segmentos_tricep_5lb(i,2));
    rms_tricep_5lb(i) = rms_valor(segmento_tricep);
%     fprintf('RMS Tríceps Segmento %d: %f\n', i, rms_tricep_5lb(i));
end

% Promedio y desviación estándar para 5 lb
mean_rms_bicep_5lb = mean(rms_bicep_5lb);
std_rms_bicep_5lb = std(rms_bicep_5lb);
mean_rms_tricep_5lb = mean(rms_tricep_5lb);
std_rms_tricep_5lb = std(rms_tricep_5lb);

fprintf('Promedio RMS Bíceps 5lb: %f, Desviación Estándar: %f\n', mean_rms_bicep_5lb, std_rms_bicep_5lb);
fprintf('Promedio RMS Tríceps 5lb: %f, Desviación Estándar: %f\n', mean_rms_tricep_5lb, std_rms_tricep_5lb);

% Segmentos para 13 lb
disp('RMS para 13 lb:');
for i = 1:5
    % Bíceps
    segmento_bicep = bicep_13lb(segmentos_bicep_13lb(i,1):segmentos_bicep_13lb(i,2));
    rms_bicep_13lb(i) = rms_valor(segmento_bicep);
%     fprintf('RMS Bíceps Segmento %d: %f\n', i, rms_bicep_13lb(i));
    
    % Tríceps
    segmento_tricep = tricep_13lb(segmentos_tricep_13lb(i,1):segmentos_tricep_13lb(i,2));
    rms_tricep_13lb(i) = rms_valor(segmento_tricep);
%     fprintf('RMS Tríceps Segmento %d: %f\n', i, rms_tricep_13lb(i));
end

% Promedio y desviación estándar para 13 lb
mean_rms_bicep_13lb = mean(rms_bicep_13lb);
std_rms_bicep_13lb = std(rms_bicep_13lb);
mean_rms_tricep_13lb = mean(rms_tricep_13lb);
std_rms_tricep_13lb = std(rms_tricep_13lb);

fprintf('Promedio RMS Bíceps 13lb: %f, Desviación Estándar: %f\n', mean_rms_bicep_13lb, std_rms_bicep_13lb);
fprintf('Promedio RMS Tríceps 13lb: %f, Desviación Estándar: %f\n', mean_rms_tricep_13lb, std_rms_tricep_13lb);

% Pesos en el eje x
pesos = [3, 5, 13];

% Promedios de RMS y desviaciones estándar ya calculados
rms_prom_bicep = [mean_rms_bicep_3lb, mean_rms_bicep_5lb, mean_rms_bicep_13lb];
std_rms_bicep = [std_rms_bicep_3lb, std_rms_bicep_5lb, std_rms_bicep_13lb];

rms_prom_tricep = [mean_rms_tricep_3lb, mean_rms_tricep_5lb, mean_rms_tricep_13lb];
std_rms_tricep = [std_rms_tricep_3lb, std_rms_tricep_5lb, std_rms_tricep_13lb];

% Subgráfica para RMS del bíceps
figure;
subplot(2, 1, 1); 
hold on;
fill([pesos, fliplr(pesos)], [rms_prom_bicep + std_rms_bicep, fliplr(rms_prom_bicep - std_rms_bicep)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(pesos, rms_prom_bicep, '-o', 'LineWidth', 2, 'Color', 'b');
title('RMS del Bíceps');
xlabel('Peso [lb]');
ylabel('RMS');
grid on;
hold off;

% Subgráfica para RMS del tríceps
subplot(2, 1, 2); 
hold on;
fill([pesos, fliplr(pesos)], [rms_prom_tricep + std_rms_tricep, fliplr(rms_prom_tricep - std_rms_tricep)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(pesos, rms_prom_tricep, '-o', 'LineWidth', 2, 'Color', 'r');
title('RMS del Tríceps');
xlabel('Peso [lb]');
ylabel('RMS');
grid on;
hold off;

sgtitle('Comparación de RMS con Desviación Estándar para Bíceps y Tríceps');

%% CALCULOS DE AUC Y STD

% Calcular el AUC para los segmentos de 3 lb
for i = 1:5
    auc_bicep_3lb(i) = trapz(t_3lb(segmentos_bicep_3lb(i,1):segmentos_bicep_3lb(i,2)), ...
                             abs(bicep_3lb(segmentos_bicep_3lb(i,1):segmentos_bicep_3lb(i,2))));
    auc_tricep_3lb(i) = trapz(t_3lb(segmentos_tricep_3lb(i,1):segmentos_tricep_3lb(i,2)), ...
                              abs(tricep_3lb(segmentos_tricep_3lb(i,1):segmentos_tricep_3lb(i,2))));
end

% Calcular el AUC para los segmentos de 5 lb
for i = 1:5
    auc_bicep_5lb(i) = trapz(t_5lb(segmentos_bicep_5lb(i,1):segmentos_bicep_5lb(i,2)), ...
                             abs(bicep_5lb(segmentos_bicep_5lb(i,1):segmentos_bicep_5lb(i,2))));
    auc_tricep_5lb(i) = trapz(t_5lb(segmentos_tricep_5lb(i,1):segmentos_tricep_5lb(i,2)), ...
                              abs(tricep_5lb(segmentos_tricep_5lb(i,1):segmentos_tricep_5lb(i,2))));
end

% Calcular el AUC para los segmentos de 13 lb
for i = 1:5
    auc_bicep_13lb(i) = trapz(t_13lb(segmentos_bicep_13lb(i,1):segmentos_bicep_13lb(i,2)), ...
                              abs(bicep_13lb(segmentos_bicep_13lb(i,1):segmentos_bicep_13lb(i,2))));
    auc_tricep_13lb(i) = trapz(t_13lb(segmentos_tricep_13lb(i,1):segmentos_tricep_13lb(i,2)), ...
                               abs(tricep_13lb(segmentos_tricep_13lb(i,1):segmentos_tricep_13lb(i,2))));
end

%STD para bicep de auc para los 3 pesos
std_auc_bicep_3lb = std(auc_bicep_3lb);
std_auc_bicep_5lb = std(auc_bicep_5lb);
std_auc_bicep_13lb = std(auc_bicep_13lb);

%STD para tricep de auc para los 3 pesos
std_auc_tricep_3lb = std(auc_tricep_3lb);
std_auc_tricep_5lb = std(auc_tricep_5lb);
std_auc_tricep_13lb = std(auc_tricep_13lb);

% Promediar el AUC de los 5 segmentos para cada peso
prom_auc_bicep_3lb = mean(auc_bicep_3lb);
prom_auc_tricep_3lb = mean(auc_tricep_3lb);

prom_auc_bicep_5lb = mean(auc_bicep_5lb);
prom_auc_tricep_5lb = mean(auc_tricep_5lb);

prom_auc_bicep_13lb = mean(auc_bicep_13lb);
prom_auc_tricep_13lb = mean(auc_tricep_13lb);

% Mostrar los resultados
disp('Área bajo la curva (AUC) para Bíceps (3 lb, 5 lb, 13 lb):');
disp(['Promedio AUC 3 lb: ', num2str(prom_auc_bicep_3lb)]);
disp(['Promedio AUC 5 lb: ', num2str(prom_auc_bicep_5lb)]);
disp(['Promedio AUC 13 lb: ', num2str(prom_auc_bicep_13lb)]);

disp('Área bajo la curva (AUC) para Tríceps (3 lb, 5 lb, 13 lb):');
disp(['Promedio AUC 3 lb: ', num2str(prom_auc_tricep_3lb)]);
disp(['Promedio AUC 5 lb: ', num2str(prom_auc_tricep_5lb)]);
disp(['Promedio AUC 13 lb: ', num2str(prom_auc_tricep_13lb)]);

% Vector con los promedios de AUC para bíceps y tríceps
prom_auc_bicep = [prom_auc_bicep_3lb, prom_auc_bicep_5lb, prom_auc_bicep_13lb];
prom_auc_tricep = [prom_auc_tricep_3lb, prom_auc_tricep_5lb, prom_auc_tricep_13lb];

% Desviaciones estándar de AUC para bíceps y tríceps
std_bicep = [std_auc_bicep_3lb, std_auc_bicep_5lb, std_auc_bicep_13lb];
std_tricep = [std_auc_tricep_3lb, std_auc_tricep_5lb, std_auc_tricep_13lb];

% Subgráfica para AUC de bíceps
figure;
subplot(2, 1, 1);
hold on;
fill([pesos, fliplr(pesos)], [prom_auc_bicep + std_bicep, fliplr(prom_auc_bicep - std_bicep)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');  
plot(pesos, prom_auc_bicep, '-o', 'LineWidth', 2, 'Color', 'b');  % Bolitas en los promedios
title('AUC Bíceps con Desviación Estándar');
xlabel('Peso (lb)');
ylabel('AUC');
grid on;
hold off;

% Subgráfica para AUC de tríceps
subplot(2, 1, 2);
hold on;
fill([pesos, fliplr(pesos)], [prom_auc_tricep + std_tricep, fliplr(prom_auc_tricep - std_tricep)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
plot(pesos, prom_auc_tricep, '-o', 'LineWidth', 2, 'Color', 'r');  % Bolitas en los promedios
title('AUC Tríceps con Desviación Estándar');
xlabel('Peso (lb)');
ylabel('AUC');
grid on;
hold off;

sgtitle('Comparación de AUC con Desviación Estándar para Bíceps y Tríceps');

% Crear una tabla solo con RMS y std
resultados_rms_std = table(...
    {'3 lb'; '5 lb'; '13 lb'}, ... % Peso
    [rms_bicep_3lb; rms_bicep_5lb; rms_bicep_13lb], ... % RMS Bíceps
    [rms_tricep_3lb; rms_tricep_5lb; rms_tricep_13lb], ... % RMS Tríceps
    [std_rms_bicep_3lb; std_rms_bicep_5lb; std_rms_bicep_13lb], ... % STD Bíceps
    [std_rms_tricep_3lb; std_rms_tricep_5lb; std_rms_tricep_13lb], ... % STD Tríceps
    'VariableNames', {'Peso', 'RMS_Biceps', 'RMS_Triceps', 'STD_Biceps', 'STD_Triceps'} ...
);
disp(resultados_rms_std);

% Crear una tabla solo con AUC y STD
resultados_auc_std = table(...
    {'3 lb'; '5 lb'; '13 lb'}, ... % Peso
    [auc_bicep_3lb; auc_bicep_5lb; auc_bicep_13lb], ... % AUC Bíceps
    [auc_tricep_3lb; auc_tricep_5lb; auc_tricep_13lb], ... % AUC Tríceps
    [std_auc_bicep_3lb; std_auc_bicep_5lb; std_auc_bicep_13lb], ... % STD Bíceps
    [std_auc_tricep_3lb; std_auc_tricep_5lb; std_auc_tricep_13lb], ... % STD Tríceps
    'VariableNames', {'Peso', 'AUC_Biceps', 'AUC_Triceps', 'STD_Biceps', 'STD_Triceps'} ...
);
disp(resultados_auc_std);

%% FATIGA (cruces por cero)

% Definir tamaño de la ventana de 10 segundos
window_size = 10 * fs; 
num_windows = floor(length(bicep_sostenido) / window_size); 

% Inicializar matrices para guardar los resultados
promedio_cruces_bic = zeros(num_windows, 1);
promedio_cruces_tri = zeros(num_windows, 1);
tiempos = (1:num_windows)' * 10; 

% Iterar sobre cada ventana de 10 segundos (50,000 muestras)
for i = 1:num_windows
    % Definir las muestras para la ventana actual
    start_idx = (i - 1) * window_size + 1;
    end_idx = start_idx + window_size - 1;
    
    % Extraer el segmento de la señal para bíceps y tríceps
    segment_bic = bicep_sostenido(start_idx:end_idx);
    segment_tri = tricep_sostenido(start_idx:end_idx);
    
    % Calcular los cruces por cero para cada ventana
    [cruces_bic, ~] = contar_cruces_por_cero(segment_bic);
    [cruces_tri, ~] = contar_cruces_por_cero(segment_tri);
    
    % Calcular el promedio de cruces por cero en la ventana de 50 segundos
    promedio_cruces_bic(i) = cruces_bic / window_size; % Promedio por muestra
    promedio_cruces_tri(i) = cruces_tri / window_size; % Promedio por muestra
end

% Graficar el promedio de cruces por cero en función del tiempo para bíceps y tríceps
figure;

% Bíceps
subplot(2, 1, 1);
plot(tiempos, promedio_cruces_bic, '-s', 'Color', [0.1, 0.6, 0.8], 'MarkerSize', 6, 'LineWidth', 1.5); % Cambiado a un estilo de línea diferente
title('Promedio de Cruces por Cero de Bíceps por Ventana de 10s');
xlabel('Tiempo (segundos)');
ylabel('Promedio de Cruces por Cero por Muestra');
grid on;

% Tríceps
subplot(2, 1, 2);
plot(tiempos, promedio_cruces_tri, '-d', 'Color', [0.9, 0.3, 0.3], 'MarkerSize', 6, 'LineWidth', 1.5); % Cambiado a un estilo de línea diferente
title('Promedio de Cruces por Cero de Tríceps por Ventana de 10s');
xlabel('Tiempo (segundos)');
ylabel('Promedio de Cruces por Cero por Muestra');
grid on;

sgtitle('Promedio de Cruces por Cero en Ventanas de 10s');

% Función para contar cruces por cero y obtener las posiciones
function [num_cruces, posiciones] = contar_cruces_por_cero(signal)
    cruces = diff(sign(signal)); 
    posiciones = find(cruces ~= 0); 
    num_cruces = length(posiciones); 
end
