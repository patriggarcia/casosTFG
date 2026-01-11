clc
clear all
close all
load parsed_trainIII_matlab

Lmodel = Lmodel_ini; % Initialization
Lmodel.update = 2; % Change this to 1 for EWMA and 2 for Iterative 
Lmodel.type = 1; % Change this to 1 for PCA and 2 for PLS
%Lmodel.lv = 3; % Initial number of LVs
Lmodel.prep = 1; % X-block prepr. 0: None, 1: Mean-center, 2: Auto-scaling 
%Lmodel.prepy = 2; % Y-block prepr. 0: None, 1: Mean-center, 2: Auto-scaling
Lmodel.nc = 50; % Number of clusters
%Lmodel.var_l = var_l; %nombre de las variables está en var_l
Lmodel.lvs = 1:2;
X = x;
obs_l = obs_l';

lambda = 1-1e-4; % Forgetting factor in EWMA
step = 1;

%% Model building (EWMA or Iterative) Iterative in this case



v=unique([0:100:size(x,1) size(x,1)]); %%% X por x
list = {};
for i=1:length(v)-1,
    list(i).x = x(v(i)+1:v(i+1),:);  %%% X por x
    list(i).obs_l = obs_l(v(i)+1:v(i+1));
end

if Lmodel.update == 1
    Lmodel = update_ewma(list, '', Lmodel, lambda, step, 1);
else
     Lmodel = update_iterative(list,'',Lmodel,step,1,1); % Iterative
end



%% Data Analysis

if Lmodel.type==2, % for PLS
    
    % Score plot
    scores_Lpls(Lmodel);   
        
    % MEDA
    map = meda_Lpls(Lmodel,0.1,111); 
    
    % reorder variables
    [map,ind] = seriation(map);
    Lmodel.XX = Lmodel.XX(ind,ind);
    Lmodel.XY = Lmodel.XY(ind,:);
    Lmodel.centr = Lmodel.centr(:,ind);
    Lmodel.var_l = Lmodel.var_l(ind);

    % oMEDAs
    dummy = zeros(100,1); % Comparison between classes 1 and 19
    dummy(find(Lmodel.class==1))=1;
    dummy(find(Lmodel.class==19))=-1;
    omeda_Lpls(Lmodel,Lmodel.centr,dummy,1);

    dummy = zeros(100,1); % Comparison between classes 1 and 11
    dummy(find(Lmodel.class==1))=1;
    dummy(find(Lmodel.class==11))=-1;
    omeda_Lpls(Lmodel,Lmodel.centr,dummy,1);
    
else %for PCA
    %% 
    load Lmodel_runtrain(4).mat
    
    % Score plot
    scores_Lpca(Lmodel);

    obs_l = obs_l';

    [Dst,Qst,Dstt,Qstt,UCLd,UCLq] = mspc_Lpca(Lmodel); %%% he añadido [lo de antes del igual, y el igual)
    
    n_PCs = 1; %cambiado porque en la gráfica el primer mínimo es en 1
    Lmodel.lvs = 1:3;

    %mspc_Lpca(Lmodel);
    x_var = var_Lpca(Lmodel);
    tscore = ((1 - x_var(n_PCs)) * Dst)/UCLd(1) + (x_var(n_PCs) * Qst)/UCLq(1);  % Compute tscore %%lo añado

    [risk, triage] = sort(tscore, 'descend');   % Triage observation based on their tscore value   %%lo añado


%%%%%%%%%%%%%%%%% voy a probar este código para poder obtener las variables
%%%%%%%%%%%%%%%%% más anómalas de las anomalias, es decir, de las
%%%%%%%%%%%%%%%%% observaciones más anómalas las variables que hacen que lo
%%%%%%%%%%%%%%%%% sean



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % MEDA
    map = meda_Lpca(Lmodel,0.1,111);
    
    % reorder variables
    [map,ind] = seriation(map);
    Lmodel.XX = Lmodel.XX(ind,ind);
    Lmodel.centr = Lmodel.centr(:,ind);
    Lmodel.centr = Lmodel.centr(:,ind);
    Lmodel.var_l = var_l;
    Lmodel.var_l = Lmodel.var_l(ind);
    
    % oMEDAs
    dummy = zeros(100,1); % Comparison between classes 1 and 19
    dummy(find(Lmodel.class==1))=1;
    dummy(find(Lmodel.class==19))=-1;
   %omeda_Lpca(Lmodel,Lmodel.centr,dummy,1);

    dummy = zeros(100,1); % Comparison between classes 1 and 11
    dummy(find(Lmodel.class==1))=1;
    dummy(find(Lmodel.class==11))=-1;
    %%omeda_Lpca(Lmodel,Lmodel.centr,dummy,1);

end

%unavez vista en el gráfic d-st vs Q-st, vamos a obtener las features más
%significativas de los outliers

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