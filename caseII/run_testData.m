% codified by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 13/Oct/16.

%% Inicia

clc
clear
load UGR16v1

Lmodel = Lmodel_ini; % Initialization
Lmodel.update = 2; % Change this to 1 for EWMA and 2 for Iterative
Lmodel.type = 1; % Change this to 1 for PCA and 2 for PLS
Lmodel.lvs = 1:3; % Max number of LVs
Lmodel.prep = 1; % X-block prepr. 0: None, 1: Mean-center, 2: Auto-scaling 
Lmodel.nc = 100; % Number of clusters
%Lmodel.path = s2wayMCpath;

step = 1;

v=unique([0:100:size(X,1) size(X,1)]);
list = {};
for i=1:length(v)-1,
    list(i).X = X(v(i)+1:v(i+1),:);
    list(i).obs_l = obs_l(v(i)+1:v(i+1));
end

Lmodel = update_iterative(list,'',Lmodel,step,1,1);

save UGR16v1 Lmodel -APPEND

%% Selección de parámetros

load UGR16v1

Lmodel.lvs = 1:20;

x_var = var_Lpca(Lmodel);

Lmodel.lvs = 1:3;

mspc_Lpca(Lmodel);

save UGR16v1 Lmodel  x_var -APPEND


%% MSNM Calibración

load UGR16v1

[LbRc,QbRc,Dstt,Qstt,UCLd,UCLq] = mspc_Lpca(Lmodel,test,0,[],[],[],[],1);

%load 'output-20160318t1047--to--output-20160705t2301';
X(:,139:end) = [];

[LbRc,QbRc,Dstc,Qstc] = mspc_Lpca(Lmodel,X,0,[],[],[],[],1);

save UGR16v1 LbRc QbRc Dstc Qstc Dstt Qstt UCLd UCLq -APPEND


%% Pred

delWE = false; % delete weekends


pred2wMC_D = Dstt;
pred2wMC_Q = Qstt;
pred2wMC = 3/138*Dstt/UCLd + 135/138*Qstt/UCLq;
predC2wMC = 3/138*Dstc/UCLd + 135/138*Qstc/UCLq;
ini = find(strcmp(obs_lt,'201607280000'));
fin = find(strcmp(obs_lt,'201608262359'));
pred2wMC_D = pred2wMC_D(ini:fin);
pred2wMC_Q = pred2wMC_Q(ini:fin);
pred2wMC = pred2wMC(ini:fin);
if delWE, classWDt = classWDt(ini:end); end

%if classWDt == classWDt(ini:end); end

save UGR16v1  pred2wMC predC2wMC  -APPEND

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

%% Evaluate: ROC (2way y 3way)

load UGR16v1


[XRoc2wMC,YRoc2wMC,TRoc2wMC,AUC2wMC] = perfcurve(Y,pred2wMC,1);
[XRocsvmMC,YRocsvmMC,TRocsvmMC,AUCsvmMC] = perfcurve(Y,predsvm,1); % lineal con inversión de plano, no lo usamos por mal definido


figure, plot(XRoc2wMC(1:1000:end),YRoc2wMC(1:1000:end)); hold on,
plot(XRocsvmMC,YRocsvmMC);
legend(sprintf('MSNM2w (AUC=%.2f)',AUC2wMC),sprintf('OCSVM (AUC=%.2f)',AUCsvmMC),'Location','southeast')
xlabel('False Positive Rate','FontSize',16); ylabel('True Positive Rate','FontSize',16)
%saveas(gcf,'figuras/ROCMC','fig'); saveas(gcf,'figuras/ROCMC.eps','epsc');




%% Visualize results per type

load UGR16v1

pred = [];
pred(:,1) = pred2wMC(1:length(y1t));
pred(:,2) = predsvm(1:length(y1t));

tity = {'MSNM','','OCSVM'};


figure,
n=30;

i=1
[kk,ord] = sort(pred(:,i),'descend');
indi = ord(1:n);
indi2 = find(sum(yt(indi,[1 2 5]),2));

subplot(2,1,i), semilogy(pred(:,i)), hold on, semilogy(indi,pred(indi,i),'or'), semilogy(indi(indi2),pred(indi(indi2),i),'sk','MarkerSize',12), axis tight
ylabel(tity(i),'Fontsize',16);

i=3;
[kk,ord] = sort(pred(:,i),'descend');
indi = ord(1:n);
indi2 = find(sum(yt(indi,[1 2 5]),2));
subplot(2,1,2), plot(pred(:,i)), hold on, plot(indi,pred(indi,i),'or'), plot(indi(indi2),pred(indi(indi2),i),'sk','MarkerSize',12), axis tight
ylabel(tity(i),'Fontsize',16);

%saveas(gcf,'figuras/data31_tmc','fig'); saveas(gcf,'figuras/data31_tmc.eps','epsc');

%% Diagnosis

load UGR16v1

clc

ini = find(strcmp(obs_lt,'201607280000'));
fin = find(strcmp(obs_lt,'201608262359'));
testo = test(ini:fin,:);
obs_lto = obs_lt(ini:fin);
classMt = classMt(ini:fin);

% 20160801t0410-0414

t = find(strcmp(obs_lto,'201608010410'))
t = t:t+4;

Lmodel.lvs = 1:rank(Lmodel.XX);
dummy = zeros(size(testo,1),1);
dummy(t) = 1; 
omeda_vec = omeda_Lpca(Lmodel,testo,dummy);
[kk1,ind]=sort(abs(omeda_vec),'descend');
C2way_1 = var_l(ind(1:11))
omeda_vec(ind(1:11))'



% 20160806t2039

t = find(strcmp(obs_lto,'201608062039'))
tf = find(strcmp(obs_lto,'201608070559'))
t = t:tf;
Lmodel.lvs = 1:rank(Lmodel.XX);
dummy = zeros(size(testo,1),1);
dummy(t) = 1; 
omeda_vec = omeda_Lpca(Lmodel,testo,dummy);
[kk2,ind]=sort(abs(omeda_vec),'descend');
C2way_2 = var_l(ind(1:10))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  representación Tscore de DataTest

%% Triage anomalies

weight_usage = true; % change this not to use (expert defined) weights 

%if weight_usage % preprocess data
  xcs = preprocess2Di(X,2);    
%else
 % xcs = preprocess2Di(X,2);
%end

%% Selection of PCs
load UGR16v1

x_var = var_Lpca(Lmodel,1);    % Plot residual variance vs # of PCs
xcs = preprocess2Di(X,2);    

n_PCs = (1:4);    
n_feat = size(xcs, 2);

Dst  =LbRc;
Qst = QbRc;


tscore = ((x_var(n_PCs)') .* Dstt)/UCLd(1) + ((1 - x_var(n_PCs)') .* Qstt)/UCLq(1);  % Compute tscore

tscoreREP = (tscore(1:24000)); 

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
topN = 1;  % número de observaciones a resaltar
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

