clear all
close

load parsed_trainIII_matlab.mat
X = x;

load parsed_testAttack_matlab.mat
obs_lt = obs_l;

load Lmodel_runtrain(4).mat



[LbRc,QbRc,Dstt,Qstt,UCLd,UCLq] = mspc_Lpca(Lmodel,x,0,[],[],[],[],1); % en test están las observaciones del conjunto de test

[LbRc,QbRc,Dstc,Qstc] = mspc_Lpca(Lmodel,X,0,[],[],[],[],1);

pred2wMC_D = Dstt;
pred2wMC_Q = Qstt;
pred2wMC = 3/138*Dstt/UCLd + 135/138*Qstt/UCLq;
predC2wMC = 3/138*Dstc/UCLd + 135/138*Qstc/UCLq;
ini = find(strcmp(obs_lt,'202407080000'));
fin = find(strcmp(obs_lt,'202407192342'));
pred2wMC_D = pred2wMC_D(ini:fin);
pred2wMC_Q = pred2wMC_Q(ini:fin);
pred2wMC = pred2wMC(ini:fin);


weight_usage = true;


x_var = var_Lpca(Lmodel,1);    % Plot residual variance vs # of PCs
xcs = preprocess2Di(X,2);    

n_PCs = (1:4);    
n_feat = size(xcs, 2);

Dst  =LbRc;
Qst = QbRc;


tscore = ((x_var(n_PCs)') .* Dstt)/UCLd(1) + ((1 - x_var(n_PCs)') .* Qstt)/UCLq(1);  % Compute tscore

tscoreREP = (tscore(1:10003)); 

[risk, triage] = sort((tscore), 'descend');   % Triage observation based on their tscore value

figure; %% Representación del Tscore vs Observación (solo top 5)
plot(((tscoreREP)), 'LineWidth', 1.5, 'Color', [0.1 0.4 0.8]); % línea azul sin marcadores
ax = gca;
ax.YScale = 'log';
hold on;

xlabel('Observation');
ylabel('Tscore');
title('Tscore vs Observation');
grid on;

% Encontrar las 5 observaciones con Tscore más alto
[sortedT, sortedIdx] = sort(tscoreREP, 'descend');
topN = 4;  % número de observaciones a resaltar
topIdx = sortedIdx(1:topN);
topVal = sortedT(1:topN);

% Resaltar los puntos con un marcador y texto
plot(topIdx, topVal, 'bo', 'MarkerFaceColor', 'w', 'MarkerSize', 1);

% Añadir etiquetas con número de observación y valor
for i = 1:topN
    text(topIdx(i), topVal(i), ...
        sprintf('obs %d (%.2f)', topIdx(i), topVal(i)), ...
        'FontSize', 9, 'Color', 'k', 'VerticalAlignment', 'bottom');
end

%%%%%%%%%%%%%%%%%%%%%%%% paso pre-diagnosis

ini = find(strcmp(obs_lt,'202408010000'));
fin = find(strcmp(obs_lt,'202408072242'));
testo = x(ini:fin,:);
obs_lto = obs_lt(ini:fin);

% 202408072239

t = find(strcmp(obs_lto,'202408072239'));


Lmodel.lvs = 1:rank(Lmodel.XX);
dummy = zeros(size(testo,1),1);
dummy(t) = 1; 
omeda_vec = omeda_Lpca(Lmodel,testo,dummy);
[kk1,ind]=sort(abs(omeda_vec),'descend');
C2way_1 = var_l(ind(1:11));


t = find(strcmp(obs_lto,'202408011341'));
t = t:t+2;


Lmodel.lvs = 1:rank(Lmodel.XX);
dummy = zeros(size(testo,1),1);
dummy(t) = 1; 
omeda_vec = omeda_Lpca(Lmodel,testo,dummy);
[kk1,ind]=sort(abs(omeda_vec),'descend');
C2way_2 = var_l(ind(1:11));


%%%%%%%%%%%%%%%%% paso deparsing

deparsing_input = fopen('deparsing_input','w+');
variables = C2way_1;
timestamp = 202408072239;

% Write variables
fwrite(deparsing_input,sprintf('features:\n{\n'));
counter = 0;
i=1;
for variable = string(variables)
  counter=counter+1;

fwrite(deparsing_input, sprintf('  [%d,%d] = %s\n', i, counter, variable)); %línea modifcada

end
fwrite(deparsing_input, sprintf('}\n'));   %línea modificada


% Write timestamps
fwrite(deparsing_input, sprintf('timestamps:\n{\n'));  %línea modificada

ts = datetime(num2str(timestamp), 'InputFormat', 'yyyyMMddHHmm');

ts = datestr(ts, 'yyyy-mm-dd HH:MM:SS');%añadido


fwrite(deparsing_input, sprintf('  [%d,%d] =  %s\n', 1, 1, ts));
fwrite(deparsing_input, sprintf('}\n'));  %línea modificada


fclose(deparsing_input);

%%%%

deparsing_input = fopen('deparsing_input2','w+');
variables = C2way_1;
timestamp = 202408011341;

% Write variables
fwrite(deparsing_input,sprintf('features:\n{\n'));
counter = 0;
i=1;
for variable = string(variables)
  counter=counter+1;

fwrite(deparsing_input, sprintf('  [%d,%d] = %s\n', i, counter, variable)); %línea modifcada

end
fwrite(deparsing_input, sprintf('}\n'));   %línea modificada


% Write timestamps
fwrite(deparsing_input, sprintf('timestamps:\n{\n'));  %línea modificada

ts = datetime(num2str(timestamp), 'InputFormat', 'yyyyMMddHHmm');

ts = datestr(ts, 'yyyy-mm-dd HH:MM:SS');%añadido


fwrite(deparsing_input, sprintf('  [%d,%d] =  %s\n', 1, 1, ts));
fwrite(deparsing_input, sprintf('}\n'));  %línea modificada


fclose(deparsing_input);
