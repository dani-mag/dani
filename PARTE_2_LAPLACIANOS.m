
%% SUPINO
data1=load('supino_laplacianos.mat');
data2=load('SUPINO_BARBIEDREAM.mat');
data3=load('Seneal_Kuakos_Supino.mat');
data4=load('SabritonesAcostado.mat');
data5=load('decubito_Sebas_potroingenieros.mat');
data6=load('DECUBITO_equipo!!.mat');

tiempo=300;  % Tiempo en segundos
tiempo2=40;
fs1=500;
muestras=fs1*tiempo;
muestras2=fs1*tiempo2;
%DATOS DE CADA EQUIPO SUPINO
ECG_EQ1 = data1.data(1:muestras); 
ECG_EQ2 = data2.data(1:muestras);
ECG_EQ3 = data3.data(1:muestras); 
ECG_EQ4 = data4.data(1:muestras);
ECG_EQ5 = data5.data(1:muestras); 
ECG_EQ6 = data6.data(1:muestras); 

%tiempo
t = (0:length(ECG_EQ1)-1) / fs1;
t2= (20:(length(ECG_EQ1)-1)+20) / fs1;
%Encontrar los picos
[pks_EQ1_50, locs_EQ1_50] = findpeaks(ECG_EQ1, 'MinPeakHeight',0.1);
[pks_EQ2_50, locs_EQ2_50] = findpeaks(ECG_EQ2, 'MinPeakHeight',0.2);
[pks_EQ3_50, locs_EQ3_50] = findpeaks(ECG_EQ3, 'MinPeakHeight',0.1);
[pks_EQ4_50, locs_EQ4_50] = findpeaks(ECG_EQ4, 'MinPeakHeight',0.1);
[pks_EQ5_50, locs_EQ5_50] = findpeaks(ECG_EQ5, 'MinPeakHeight',0.1);
[pks_EQ6_50, locs_EQ6_50] = findpeaks(ECG_EQ6, 'MinPeakDistance',330);

%Obtener los intervalos RR
RR_intervals_EQ1_50 = diff(locs_EQ1_50) / fs1 * 1000;
RR_intervals_EQ2_50 = diff(locs_EQ2_50) / fs1 * 1000;
RR_intervals_EQ3_50 = diff(locs_EQ3_50) / fs1 * 1000;
RR_intervals_EQ4_50 = diff(locs_EQ4_50) / fs1 * 1000;
RR_intervals_EQ5_50 = diff(locs_EQ5_50) / fs1 * 1000;
RR_intervals_EQ6_50 = diff(locs_EQ6_50) / fs1 * 1000;

%Obtener SDNN para cada señal
SDNN_EQ1_50 = std(RR_intervals_EQ1_50);
SDNN_EQ2_50 = std(RR_intervals_EQ2_50);
SDNN_EQ3_50 = std(RR_intervals_EQ3_50);
SDNN_EQ4_50 = std(RR_intervals_EQ4_50);
SDNN_EQ5_50 = std(RR_intervals_EQ5_50);
SDNN_EQ6_50 = std(RR_intervals_EQ6_50);

% Calcular RMSSD para cada señal
RMSSD_EQ1_50 = sqrt(mean(diff(RR_intervals_EQ1_50).^2));
RMSSD_EQ2_50 = sqrt(mean(diff(RR_intervals_EQ2_50).^2));
RMSSD_EQ3_50 = sqrt(mean(diff(RR_intervals_EQ3_50).^2));
RMSSD_EQ4_50 = sqrt(mean(diff(RR_intervals_EQ4_50).^2));
RMSSD_EQ5_50 = sqrt(mean(diff(RR_intervals_EQ5_50).^2));
RMSSD_EQ6_50 = sqrt(mean(diff(RR_intervals_EQ6_50).^2));

%MOSTRAR RESULTADOS
fprintf('SDNN para ECG supino EQ1 (1): %.2f ms\n', SDNN_EQ1_50);
fprintf('SDNN para ECG supino EQ2 (1): %.2f ms\n', SDNN_EQ2_50);
fprintf('SDNN para ECG supino EQ3 (1): %.2f ms\n', SDNN_EQ3_50);
fprintf('SDNN para ECG supino EQ4 (1): %.2f ms\n', SDNN_EQ4_50);
fprintf('SDNN para ECG supino EQ5 (1): %.2f ms\n', SDNN_EQ5_50);
fprintf('SDNN para ECG supino EQ6 (1): %.2f ms\n', SDNN_EQ6_50);

fprintf('RMSSD para ECG supino EQ1 (1): %.2f ms\n', RMSSD_EQ1_50);
fprintf('RMSSD para ECG supino EQ2 (1): %.2f ms\n', RMSSD_EQ2_50);
fprintf('RMSSD para ECG supino EQ3 (1): %.2f ms\n', RMSSD_EQ3_50);
fprintf('RMSSD para ECG supino EQ4 (1): %.2f ms\n', RMSSD_EQ4_50);
fprintf('RMSSD para ECG supino EQ5 (1): %.2f ms\n', RMSSD_EQ5_50);
fprintf('RMSSD para ECG supino EQ6 (1): %.2f ms\n', RMSSD_EQ6_50);

%Plotear
% Plotear y marcar los picos R con una 'x' en cada señal
figure;
plot(t, ECG_EQ1);
hold on;
plot(t(locs_EQ1_50), ECG_EQ1(locs_EQ1_50), 'x', 'MarkerSize', 10, 'Color', 'r');  % Marcar picos R con 'x'
title('ECG EQ1');
hold off;

figure;
plot(t, ECG_EQ2);
hold on;
plot(t(locs_EQ2_50), ECG_EQ2(locs_EQ2_50), 'x', 'MarkerSize', 10, 'Color', 'r');  % Marcar picos R con 'x'
title('ECG EQ2');
hold off;

figure;
plot(t, ECG_EQ3);
hold on;
plot(t(locs_EQ3_50), ECG_EQ3(locs_EQ3_50), 'x', 'MarkerSize', 10, 'Color', 'r');  % Marcar picos R con 'x'
title('ECG EQ3');
hold off;

figure;
plot(t, ECG_EQ4);
hold on;
plot(t(locs_EQ4_50), ECG_EQ4(locs_EQ4_50), 'x', 'MarkerSize', 10, 'Color', 'r');  % Marcar picos R con 'x'
title('ECG EQ4');
hold off;

figure;
plot(t, ECG_EQ5);
hold on;
plot(t(locs_EQ5_50), ECG_EQ5(locs_EQ5_50), 'x', 'MarkerSize', 10, 'Color', 'r');  % Marcar picos R con 'x'
title('ECG EQ5');
hold off;

figure;
plot(t2, ECG_EQ6);
hold on;
plot(t2(locs_EQ6_50), ECG_EQ6(locs_EQ6_50), 'x', 'MarkerSize', 10, 'Color', 'r');  % Marcar picos R con 'x'
title('ECG EQ6');
hold off;

% Organizar los valores de SDNN en una matriz 
SDNN_values = [SDNN_EQ1_50; SDNN_EQ2_50; SDNN_EQ3_50; SDNN_EQ4_50; SDNN_EQ5_50; SDNN_EQ6_50];

RMSSD_values= [RMSSD_EQ1_50; RMSSD_EQ2_50; RMSSD_EQ3_50; RMSSD_EQ4_50; RMSSD_EQ5_50; RMSSD_EQ6_50];

% Crear el boxplot para comparar los SDNNN acostados
figure;
boxplot(SDNN_values);
xlabel('Equipos');
ylabel('SDNN (ms)');
title('Comparación de SDNN supino entre Equipos');
xticklabels({'Equipo 1', 'Equipo 2', 'Equipo 3', 'Equipo 4', 'Equipo 5', 'Equipo 6'});

% Crear el boxplot para comparar los RMSSD acostados
figure;
boxplot(RMSSD_values);
xlabel('Equipos');
ylabel('RMSSD (ms)');
title('Comparación de RMSSD supino entre Equipos');
xticklabels({'Equipo 1', 'Equipo 2', 'Equipo 3', 'Equipo 4', 'Equipo 5', 'Equipo 6'});

%% ORTOSTÁTICOS
data1_2=load('pie_laplacianos.mat');
data2_2=load('PARADA_BARBIEDREAM.mat');
data3_2=load('Seneal_Kuakos_Pie.mat');
data4_2=load('SabritonesParado.mat');
data5_2=load('parado_Sebas_potroingenieros.mat');
data6_2=load('PARADO_equipo!!.mat');

%DATOS DE CADA EQUIPO PARADOS
ECG_EQ1P = data1_2.data(1:muestras); 
ECG_EQ2P = data2_2.data(1:muestras);
ECG_EQ3P = data3_2.data(1:muestras); 
ECG_EQ4P = data4_2.data(1:muestras);
ECG_EQ5P = data5_2.data(1:muestras); 
ECG_EQ6P = data6_2.data(1:muestras); 

%Encontrar los picos
[pks_EQ1P_50, locs_EQ1P_50] = findpeaks(ECG_EQ1P, 'MinPeakHeight',0.1);
[pks_EQ2P_50, locs_EQ2P_50] = findpeaks(ECG_EQ2P, 'MinPeakHeight',0.2);
[pks_EQ3P_50, locs_EQ3P_50] = findpeaks(ECG_EQ3P, 'MinPeakDistance',250);
[pks_EQ4P_50, locs_EQ4P_50] = findpeaks(ECG_EQ4P, 'MinPeakDistance',272);
[pks_EQ5P_50, locs_EQ5P_50] = findpeaks(ECG_EQ5P, 'MinPeakDistance',217);
[pks_EQ6P_50, locs_EQ6P_50] = findpeaks(ECG_EQ6P, 'MinPeakDistance',207.2);

%Obtener los intervalos RR
RR_intervals_EQ1P_50 = diff(locs_EQ1P_50) / fs1 * 1000;
RR_intervals_EQ2P_50 = diff(locs_EQ2P_50) / fs1 * 1000;
RR_intervals_EQ3P_50 = diff(locs_EQ3P_50) / fs1 * 1000;
RR_intervals_EQ4P_50 = diff(locs_EQ4P_50) / fs1 * 1000;
RR_intervals_EQ5P_50 = diff(locs_EQ5P_50) / fs1 * 1000;
RR_intervals_EQ6P_50 = diff(locs_EQ6P_50) / fs1 * 1000;

%Obtener SDNN para cada señal
SDNN_EQ1P_50 = std(RR_intervals_EQ1P_50);
SDNN_EQ2P_50 = std(RR_intervals_EQ2P_50);
SDNN_EQ3P_50 = std(RR_intervals_EQ3P_50);
SDNN_EQ4P_50 = std(RR_intervals_EQ4P_50);
SDNN_EQ5P_50 = std(RR_intervals_EQ5P_50);
SDNN_EQ6P_50 = std(RR_intervals_EQ6P_50);

% Calcular RMSSD para cada señal
RMSSD_EQ1P_50 = sqrt(mean(diff(RR_intervals_EQ1P_50).^2));
RMSSD_EQ2P_50 = sqrt(mean(diff(RR_intervals_EQ2P_50).^2));
RMSSD_EQ3P_50 = sqrt(mean(diff(RR_intervals_EQ3P_50).^2));
RMSSD_EQ4P_50 = sqrt(mean(diff(RR_intervals_EQ4P_50).^2));
RMSSD_EQ5P_50 = sqrt(mean(diff(RR_intervals_EQ5P_50).^2));
RMSSD_EQ6P_50 = sqrt(mean(diff(RR_intervals_EQ6P_50).^2));

%MOSTRAR RESULTADOS
fprintf('SDNN para ECG ortostático EQ1 (1): %.2f ms\n', SDNN_EQ1P_50);
fprintf('SDNN para ECG ortostático EQ2 (1): %.2f ms\n', SDNN_EQ2P_50);
fprintf('SDNN para ECG ortostático EQ3 (1): %.2f ms\n', SDNN_EQ3P_50);
fprintf('SDNN para ECG ortostático EQ4 (1): %.2f ms\n', SDNN_EQ4P_50);
fprintf('SDNN para ECG ortostático EQ5 (1): %.2f ms\n', SDNN_EQ5P_50);
fprintf('SDNN para ECG ortostático EQ6 (1): %.2f ms\n', SDNN_EQ6P_50);

fprintf('RMSSD para ECG ortostático EQ1 (1): %.2f ms\n', RMSSD_EQ1P_50);
fprintf('RMSSD para ECG ortostático EQ2 (1): %.2f ms\n', RMSSD_EQ2P_50);
fprintf('RMSSD para ECG ortostático EQ3 (1): %.2f ms\n', RMSSD_EQ3P_50);
fprintf('RMSSD para ECG ortostático EQ4 (1): %.2f ms\n', RMSSD_EQ4P_50);
fprintf('RMSSD para ECG ortostático EQ5 (1): %.2f ms\n', RMSSD_EQ5P_50);
fprintf('RMSSD para ECG ortostático EQ6 (1): %.2f ms\n', RMSSD_EQ6P_50);

%Plotear
% Plotear y marcar los picos R con una 'x' en cada señal
figure;
plot(t, ECG_EQ1P);
hold on;
plot(t(locs_EQ1P_50), ECG_EQ1P(locs_EQ1P_50), 'x', 'MarkerSize', 10, 'Color', 'r');  % Marcar picos R con 'x'
title('ECG EQ1');
hold off;

figure;
plot(t, ECG_EQ2P);
hold on;
plot(t(locs_EQ2P_50), ECG_EQ2P(locs_EQ2P_50), 'x', 'MarkerSize', 10, 'Color', 'r');  % Marcar picos R con 'x'
title('ECG EQ2');
hold off;

figure;
plot(t, ECG_EQ3P);
hold on;
plot(t(locs_EQ3P_50), ECG_EQ3P(locs_EQ3P_50), 'x', 'MarkerSize', 10, 'Color', 'r');  % Marcar picos R con 'x'
title('ECG EQ3');
hold off;

figure;
plot(t, ECG_EQ4P);
hold on;
plot(t(locs_EQ4P_50), ECG_EQ4P(locs_EQ4P_50), 'x', 'MarkerSize', 10, 'Color', 'r');  % Marcar picos R con 'x'
title('ECG EQ4');
hold off;

figure;
plot(t, ECG_EQ5P);
hold on;
plot(t(locs_EQ5P_50), ECG_EQ5P(locs_EQ5P_50), 'x', 'MarkerSize', 10, 'Color', 'r');  % Marcar picos R con 'x'
title('ECG EQ5');
hold off;

figure;
plot(t2, ECG_EQ6P);
hold on;
plot(t2(locs_EQ6P_50), ECG_EQ6P(locs_EQ6P_50), 'x', 'MarkerSize', 10, 'Color', 'r');  % Marcar picos R con 'x'
title('ECG EQ6');
hold off;

% Organizar los valores de SDNN en una matriz 
SDNN_values_Stand = [SDNN_EQ1P_50; SDNN_EQ2P_50; SDNN_EQ3P_50; SDNN_EQ4P_50; SDNN_EQ5P_50; SDNN_EQ6P_50];
RMSSD_values_Stand= [RMSSD_EQ1P_50; RMSSD_EQ2P_50; RMSSD_EQ3P_50; RMSSD_EQ4P_50; RMSSD_EQ5P_50; RMSSD_EQ6P_50];

% Crear el boxplot para comparar los SDNN parados
figure;
boxplot(SDNN_values_Stand);
xlabel('Equipos');
ylabel('SDNN (ms)');
title('Comparación de SDNN supino entre Equipos');
xticklabels({'Equipo 1', 'Equipo 2', 'Equipo 3', 'Equipo 4', 'Equipo 5', 'Equipo 6'});

% Crear el boxplot para comparar los RMSSD parados
figure;
boxplot(RMSSD_values_Stand);
xlabel('Equipos');
ylabel('RMSSD (ms)');
title('Comparación de RMSSD ortostático entre Equipos');
xticklabels({'Equipo 1', 'Equipo 2', 'Equipo 3', 'Equipo 4', 'Equipo 5', 'Equipo 6'});

% Organizar los valores de SDNN en una matriz (acostados y parados juntos)
SDNN_values_combined = [SDNN_values; SDNN_values_Stand];

% Crear etiquetas para las condiciones (Acostado y Parado)
group_labels = [repmat({'Supino'}, 6, 1); repmat({'Ortostático'}, 6, 1)];

% Crear el boxplot comparando SDNN para ambas posiciones en un solo gráfico
figure;
boxplot(SDNN_values_combined, group_labels);
xlabel('Posición');
ylabel('SDNN (ms)');
title('Comparación de SDNN en posición supino y ortostático');

% Organizar los valores de RMSSD en una matriz (acostados y parados juntos)
RMSSD_values_combined = [RMSSD_values; RMSSD_values_Stand];

% Crear etiquetas para las condiciones (Acostado y Parado)
group_labels_R = [repmat({'Supino'}, 6, 1); repmat({'Ortostático'}, 6, 1)];

% Crear el boxplot comparando RMSSD para ambas posiciones en un solo gráfico
figure;
boxplot(RMSSD_values_combined, group_labels_R);
xlabel('Posición');
ylabel('RMSSD (ms)');
title('Comparación de RMSSD en posición supino y ortostático');