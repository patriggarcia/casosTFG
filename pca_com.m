% PCA based Multivariate Statistical Network Monitoring
% VAST2012 Challenge dataset
% 
%  Data set and Analysis: 
% 
%   Camacho, J., Perez-Villegas, A., Garcia-Teodoro, P., Macia-Fernandez, G. 
%   PCA-based Multivariate Statistical Network Monitoring for Anomaly Detection.
%   Computers & Security, 2016, 59: 118-137. 
%
%   Camacho, J., Macia-Fernandez, G., Diaz-Verdejo, J., Garcia-Teodoro, P. 
%   Tackling the Big Data 4 Vs for Anomaly Detection. INFOCOM'2014 Workshop 
%   on Security and Privacy in Big Data, Toronto (Canada), 2014.
%
% The VAST 2012 2nd mini challenge is a benchmark for visualization in cybersecurity 
% (http://www.vacommunity.org/VAST+Challenge+2012)
% 
% The goal is to identify cybersecurity issues in the data collected during two days from a
% computer network. During those days, a number of non-legitimate programs were found
% to be running on several computers, slowing them down. A cyber-forensics operation is
% required to discover the root causes for this strange behavior.
% 
% Two typical sources of data are collected from the network: firewall and Intrusion
% Detection System (IDS) logs. The firewall analyses the incoming and outgoing data traffic
% in the network, and records in a log file all connection attempts that are blocked according
% to security policies. The IDS employs higher level intelligence to identify cybersecurity
% incidents in data traffic. It also stores the results in a log file, though it does not block any
% traffic connection. Also, typically, it only analyses a sub-set (sample) of the total traffic.
% 
% A total of 2345 observations, each one with the information for one minute, are obtained.
% For every sampling period of one minute, we defined a set of 112 variables that represent
% the information from the two data sources: 69 variables for the firewall log and 43 for
% the IDS log. The number of variables was reduced to 95 by discarding those with constant 
% value throughout the capture period. The definition of the variables is
% introduced in Tackling the Big Data 4 Vs for Anomaly Detection. INFOCOM'2014 Workshop on 
% Security and Privacy in Big Data, Toronto (Canada), 2014.

% coded by: Jose Manuel Garcia Gimenez
%           Jose Camacho Paez.
% last modification: 15/Feb/18.

%% Inicilization

clear all
close all
clc

load parsed_data_matlab.mat

weight_usage = true; % change this not to use (expert defined) weights 

if weight_usage % preprocess data
  xcs = preprocess2D(x,2,weight');    
else
  xcs = preprocess2D(x,2);
end


%% Selection of PCs

x_var = var_pca(xcs,0:20,0);    % Plot residual variance vs # of PCs

n_PCs = 9    % 9 PCs are choosen (first minimum in ckf)
n_feat = size(xcs, 2);

%% Compute D-st and Q-st in PCA-based MSPC and issue visualization

[Dst,Qst,Dstt,Qstt,UCLd,UCLq] = mspc_pca(xcs,1:n_PCs,[],0,11,[],[],[],[],1);


%% Triage anomalies


tscore = ((1 - x_var(n_PCs)) * Dst)/UCLd(1) + (x_var(n_PCs) * Qst)/UCLq(1);  % Compute tscore


[risk, triage] = sort(tscore, 'descend');   % Triage observation based on their tscore value

%% Representación del Tscore vs Observación (solo top 5)
figure;
plot(tscore, 'LineWidth', 1.5, 'Color', [0.1 0.4 0.8]); % línea azul sin marcadores
hold on;

xlabel('Observation');
ylabel('Tscore');
title('Tscore vs Observation');
grid on;

% Encontrar las 5 observaciones con Tscore más alto
[sortedT, sortedIdx] = sort(tscore, 'descend');
topN = 5;  % número de observaciones a resaltar
topIdx = sortedIdx(1:topN);
topVal = sortedT(1:topN);

% Resaltar los puntos con un marcador y texto
plot(topIdx, topVal, 'bo', 'MarkerFaceColor', 'w', 'MarkerSize', 5);

% Añadir etiquetas con número de observación y valor
for i = 1:topN
    text(topIdx(i) + 0.5, topVal(i), ...
        sprintf('obs %d (%.2f)', topIdx(i), topVal(i)), ...
        'FontSize', 9, 'Color', 'k', 'VerticalAlignment', 'bottom');
end

% (Opcional) Línea del umbral visual, si la usas en tu triage
%yline(0.25 * max(tscore), '--k', '0.25·max threshold', 'LabelHorizontalAlignment','left');

%legend({'Top 5 anomalies','Umbral'}, 'Location', 'best');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

anom = find(risk > 0.25 * max(risk));    % Select the most anomalous observation

for i=1:length(anom)
  omeda_vec = omeda_pca(xcs,1:size(xcs,2),xcs(triage(i),:),1,0,100,var_l);   % Compute and plot oMEDA vector for the detected anomalies
  variables{i} = var_l(find(omeda_vec > 0.2*max(omeda_vec))); % register variables and timestamp for de-parsing
  timestamps{i} = obs_l(triage(i)); 
  
  %%%%%%%%%%%%%%%%%%%%%%%%
% --- robust plotting según prefijo de variable ---
% omeda_vec ya calculado por omeda_pca
nVars = numel(omeda_vec);

% 1) Normalizar var_l a cell array de cadenas
if iscell(var_l)
    names = var_l(:);             % asegurar columna
elseif isstring(var_l)
    names = cellstr(var_l(:));
elseif ischar(var_l)
    % char matrix: cada fila es un nombre (o una sola cadena larga)
    names = cellstr(var_l);       % convierte filas en cellstr
else
    % fallback: intentar convertir a cellstr
    try
        names = cellstr(var_l);
    catch
        % crear nombres por defecto
        names = arrayfun(@(k) sprintf('var_%d',k), 1:nVars, 'UniformOutput', false)';
    end
end

% 2) Trim de espacios y asegurar que haya nVars nombres
names = strtrim(names);

if numel(names) > nVars
    names = names(1:nVars);   % recortar extras
elseif numel(names) < nVars
    % rellenar con nombres genericos
    addCount = nVars - numel(names);
    more = arrayfun(@(k) sprintf('var_%d', numel(names)+k), 1:addCount, 'UniformOutput', false);
    names = [names; more'];
end

% 3) Detectar prefijos (case-insensitive)
is_nf  = startsWith(names, 'fw',  'IgnoreCase', true);
is_ids = startsWith(names, 'ids', 'IgnoreCase', true);

% 4) Construir tabla de colores
col_default = [0.7 0.7 0.7];     % gris
col_nf      = [0.8500 0.3250 0.0980]; % naranja
col_ids     = [0 0.4470 0.7410]; % azul


CData = repmat(col_default, nVars, 1);
CData(is_nf, :)  = repmat(col_nf,  sum(is_nf),  1);
CData(is_ids, :) = repmat(col_ids, sum(is_ids), 1);

% 5) Dibujar barras con color por elemento
figure; 
hBar = bar(1:nVars, omeda_vec, 'FaceColor', 'flat', 'EdgeColor', 'none');
hBar.CData = CData;   % asignar colores por barra
xlim([0 nVars+1]);
xlabel('Variable index');
ylabel('oMEDA value');
title(sprintf('obs. %d  -->  Tscore: %.3g', triage(i), risk(i)));

% 6) Leyenda (si hay alguno de cada tipo)
hold on;
legendEntries = {};
legendHandles = [];
if any(is_nf)
    h = plot(nan, nan, 's', 'MarkerEdgeColor','none','MarkerFaceColor',col_nf);
    legendHandles(end+1) = h; %#ok<SAGROW>
    legendEntries{end+1} = 'fw';
end
if any(is_ids)
    h = plot(nan, nan, 's', 'MarkerEdgeColor','none','MarkerFaceColor',col_ids);
    legendHandles(end+1) = h; %#ok<SAGROW>
    legendEntries{end+1} = 'ids';
end
%{
h = plot(nan, nan, 's', 'MarkerEdgeColor','none','MarkerFaceColor',col_default);
legendHandles(end+1) = h;
legendEntries{end+1} = 'others';
%}
if ~isempty(legendHandles)
    legend(legendHandles, legendEntries, 'Location', 'best');
end


hold off;
  %%%%%%%%%%%%%%
  
  title(strcat('obs. ',num2str(triage(i)),'  --> ', ' Tscore: ', num2str(risk(i))))

  %%%%%%%
  xlabel('Variable index');
    ylabel('oMEDA value');

end


%% Store results 

%%%save "/home/mbda/VAST-MC2/anomalies_PCA" variables timestamps
save 'anomalies_PCA' variables timestamps



































