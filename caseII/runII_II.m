clear clc
clear
load caseII

Lmodel = Lmodel_ini; % Initialization
Lmodel.update = 2; % Change this to 1 for EWMA and 2 for Iterative 
Lmodel.type = 1; % Change this to 1 for PCA and 2 for PLS
Lmodel.prep = 1; % X-block prepr. 0: None, 1: Mean-center, 2: Auto-scaling 
Lmodel.nc = 100; % Number of clusters
Lmodel.lvs = 1:2;


lambda = 1-1e-4; % Forgetting factor in EWMA
step = 1;

%% Model building (EWMA or Iterative)

v=unique([0:100:size(X,1) size(X,1)]);
list = {};
for i=1:length(v)-1,
    list(i).X = X(v(i)+1:v(i+1),:);
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
    
    % Score plot
    scores_Lpca(Lmodel);

    [Dst,Qst,Dstt,Qstt,UCLd,UCLq] = mspc_Lpca(Lmodel); %%% he añadido [lo de antes del igual, y el igual)
    
    n_PCs = 4;
    Lmodel.lvs = 1:3;

    mspc_Lpca(Lmodel);
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
    omeda_Lpca(Lmodel,Lmodel.centr,dummy,1);

    dummy = zeros(100,1); % Comparison between classes 1 and 11
    dummy(find(Lmodel.class==1))=1;
    dummy(find(Lmodel.class==11))=-1;
    omeda_Lpca(Lmodel,Lmodel.centr,dummy,1);

end

%unavez vista en el gráfic d-st vs Q-st, vamos a obtener las features más
%significativas de los outliers

%% PRE-DIAGNOSIS STEP
load Lmodel_II_Train.mat
load UGR16v1.mat

ini = find(strcmp(obs_l,'201603190000'));
fin = find(strcmp(obs_l,'201606252359'));
train = X(ini:fin,:);
obs_lto = obs_l(ini:fin);
classMt = classM(ini:fin);

% 201603261229

t = find(strcmp(obs_lto,'201603261229'))
t = t:t+4;

x = X;

Lmodel.lvs = 1:rank(Lmodel.XX);
dummy = zeros(size(train,1),1);
dummy(t) = 1; 
omeda_vec = omeda_Lpca(Lmodel,train,dummy);
[kk1,ind]=sort(abs(omeda_vec),'descend');
C2way_29 = var_l(ind(1:11)); % EN ESTA VARIABLE SE ALMACENAN LAS FEATURES MÁS SIGNIFICATIVAS DEL OUTLIER, CUYO TIMESTAMP SE INTRODUCE ARRIBA
omeda_vec(ind(1:11))'


% 201603261230

t = find(strcmp(obs_lto,'201603261230'))
t = t:t+4;

x = X;

Lmodel.lvs = 1:rank(Lmodel.XX);
dummy = zeros(size(train,1),1);
dummy(t) = 1; 
omeda_vec = omeda_Lpca(Lmodel,train,dummy);
[kk2,ind]=sort(abs(omeda_vec),'descend');
C2way_30 = var_l(ind(1:11)); % EN ESTA VARIABLE SE ALMACENAN LAS FEATURES MÁS SIGNIFICATIVAS DEL OUTLIER, CUYO TIMESTAMP SE INTRODUCE ARRIBA
omeda_vec(ind(1:11))';

% 201605131003

t = find(strcmp(obs_lto,'201605131003'))
t = t:t+4;

x = X;

Lmodel.lvs = 1:rank(Lmodel.XX);
dummy = zeros(size(train,1),1);
dummy(t) = 1; 
omeda_vec = omeda_Lpca(Lmodel,train,dummy);
[kk3,ind]=sort(abs(omeda_vec),'descend');
C2way_3 = var_l(ind(1:14)); % EN ESTA VARIABLE SE ALMACENAN LAS FEATURES MÁS SIGNIFICATIVAS DEL OUTLIER, CUYO TIMESTAMP SE INTRODUCE ARRIBA
omeda_vec(ind(1:11))'



  
