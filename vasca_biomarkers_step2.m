%% Ferrer et al. INTESTINAL MICROBIOTA INFLAMMATORY PROFILE IN LATE-FETAL GROWTH RESTRICTION WOMEN, 2025. 
% Script to analyze the inflammatory biomarkers (after outliers removal)
%
% We consider the following model: 
%
%   X = A + B + C(A) + AB + E
%
% with A: Class (Control vs IUGR), B: Time (T3, Labour), C(A): Individual
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 15/May/2025
%
% Copyright (C) 2025  University of Granada, Granada
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


%% Prepare data

clear
close all
clc

% comment the line below to observe full randomness
rng(0); % we fix the random seed for repetitivity, but permutation testing is stochastic so a little variability is expected in the p-values and associated figures

pvalue = 0.05 

load biomark_vasca.mat

X = X_ASCAt;
Xraw = X;
F = F_ASCA;
var_l = var_final;

ind = find(F(:,2)==3);% do not consider newborns
X(ind,:) = [];
Xraw(ind,:) = [];
F(ind,:) = [];

indout = find(F(:,3)==16); % Deleting outlier found in first step
X(indout,:) = [];
F(indout,:) = [];

indout = find(F(:,3)==7); % Deleting outlier found in first step
X(indout,:) = [];
F(indout,:) = [];

indout = find(F(:,3)==74); % Deleting outlier found in first step
X(indout,:) = [];
F(indout,:) = [];

indout = find(F(:,3)==43); % Deleting outlier found in first step
X(indout,:) = [];
F(indout,:) = [];

indout = find(F(:,3)==31); % Deleting outlier found in first step
X(indout,:) = [];
F(indout,:) = [];

Xr =  rankTransform(X);
Xn = X; Xn(find(isnan(X)))=0;
Xrn = X; Xrn(find(isnan(X)))=0; 
X = [X Xr]; % Careful: use [X/norm(Xn) Xr/norm(Xrn)] if not autoscaled

for i=1:length(var_final), var_l{end+1} = strcat(var_l{i},'_rank'); end

clear X_ASCAt F_ASCA var_final Xr Xn Xrn


%% FDR (Q-value)
       

[T, parglmoMC] = parglmMC(X, F, 'Model', [1 2], 'Mtc', -1, 'Nested', [1 3], 'Random', [0 0 1]);
T.Source(2:5)={'F1: LFR/Ctrl','F2: Time','F3: Individual','Int Class x Time'}


%% Compare with after rank tranform (if similar, raw data preferred) 

half = length(var_l)/2;
noTr = 1:half;
Tr = half+1:2*half;

diff1 = find( (parglmoMC.p(noTr,1)<pvalue | parglmoMC.p(Tr,1)<pvalue) & abs(parglmoMC.p(noTr,1)-parglmoMC.p(Tr,1)) > pvalue);
if length(diff1)>0, disp('Check the following variables in Factor 1:'), end
disp(var_l(diff1))

diff2 = find( (parglmoMC.p(noTr,2)<pvalue | parglmoMC.p(Tr,2)<pvalue) & abs(parglmoMC.p(noTr,2)-parglmoMC.p(Tr,2)) > pvalue);
if length(diff2)>0, disp('Check the following variables in Factor 2:'), end
disp(var_l(diff2))

diff4 = find( (parglmoMC.p(noTr,4)<pvalue | parglmoMC.p(Tr,4)<pvalue) & abs(parglmoMC.p(noTr,4)-parglmoMC.p(Tr,4)) > pvalue);
if length(diff4)>0, disp('Check the following variables in the Interaction:'), end;
disp(var_l(diff4))

X2 = X(:,noTr);
var_l2 = var_l(noTr);
    
if length(diff1)==0 && length(diff2)==0 && length(diff4)==0 
    disp('Rank transformation does not make a difference, so we are safe to use the original data.'); 
else
    disp('We use the rank tranformation for the Q-value of non-normal variables.');
    parglmoMC.p(unique([diff1;diff2;diff4]),:) = parglmoMC.p(unique([diff1;diff2;diff4])+half,:);
    X2(:,unique([diff1;diff2;diff4])) = X(:,unique([diff1;diff2;diff4])+half);
    var_l2(unique([diff1;diff2;diff4])) = var_l(unique([diff1;diff2;diff4])+half);
end


%% Univariate inference: Q-value with selected features

[T, parglmoMC] = parglmMC(X2, F, 'Model', [1 2], 'Mtc', -1, 'Nested', [1 3], 'Random', [0 0 1]);
T.Source(2:5)={'F1: LFR/Ctrl','F2: Time','F3: Individual','Int LFR/Ctrl x Time'}


TQ = table(var_l2', parglmoMC.p(:,1), parglmoMC.p(:,2), parglmoMC.p(:,4),'VariableNames', {'Labels','QvalueLFRCtrl','QvalueTime','QvalueInt'})


h=figure; hold on
for factor=1:3
    plot(-log10(parglmoMC.p(parglmoMC.ordFactors(factor,:),factor)))
end
plot(-log10(parglmoMC.p(parglmoMC.ordInteractions(1,:),4)))
plot([1 size(X2,2)],-log10([pvalue pvalue]),'r--')
legend('Factor LFR/Ctrl','Factor Time','Factor Individual','LFR/Ctrl x Time','\alpha=0.05','Location','south')
a=get(h,'CurrentAxes');
set(a,'FontSize',14)
set(a,'YScale','log')
ylabel('-log_{10}(p-value)','FontSize',18)
xlabel('Variables in selected order','FontSize',18)
title('Manhattan Plot Qvalue')
axis tight

saveas(gcf,'Figures/FigManhattan_Qv_bio');
saveas(gcf,'Figures/FigManhattan_Qv_bio.jpg','jpg');
saveas(gcf,'Figures/FigManhattan_Qv_bio.png','png');


%% Multivariate inference: ASCA 

[T, parglmo] = parglm(X2, F, 'Model', [1 2], 'Nested', [1 3], 'Random', [0 0 1]);
T.Source(2:5)={'F1: LFR/Ctrl','F2: Time','F3: Individual','Int LFR/Ctrl x Time'}


%% Multivariate inference with variable selection: VASCA 

[T, parglmoVS] = parglmVS(X2, F, 'Model', [1 2], 'Nested', [1 3], 'Random', [0 0 1]);
T.Source(2:5)={'F1: LFR/Ctrl','F2: Time','F3: Individual','Int LFR/Ctrl x Time'}

TtQ = table(var_l2', parglmoVS.p(:,1), parglmoVS.p(:,2), parglmoVS.p(:,4),'VariableNames', {'Labels','VASCALFRCtrl','VASCATime','VASCAInt'})

h=figure; hold on
for factor=1:3
    plot(-log10(parglmoVS.p(parglmoVS.ordFactors(factor,:),factor)))
end
plot(-log10(parglmoVS.p(parglmoVS.ordInteractions(1,:),4)))
plot([1 size(X2,2)],-log10([pvalue pvalue]),'r--')
legend('Factor LFR/Ctrl','Factor Time','Factor Individual','LFR/Ctrl x Time','\alpha=0.05','Location','south')
a=get(h,'CurrentAxes');
set(a,'FontSize',14)
set(a,'YScale','log')
ylabel('-log_{10}(p-value)','FontSize',18)
xlabel('Variables in selected order','FontSize',18)
title('Manhattan Plot vASCA')
axis tight

saveas(gcf,'Figures/FigManhattan_bio');
saveas(gcf,'Figures/FigManhattan_bio.jpg','jpg');
saveas(gcf,'Figures/FigManhattan_bio.png','png');


%% Display results: VASCA

vascao = vasca(parglmoVS, pvalue);


% LFR/Ctrl factor

i=1;
if isfield(vascao.factors{i},'scores')

    TvASCA = table((1:length(vascao.factors{i}.ind))',var_l2(sort(vascao.factors{i}.ind))', vascao.p(sort(vascao.factors{i}.ind),i),'VariableNames', {'Order','Label',sprintf('PvalueF%i',i)})
    
    scores(vascao.factors{i}, 'Title', 'Case/control', 'ObsClass', vascao.design(:,i));
    legend('Control', 'LFG');


    if length(vascao.factors{i}.ind)==1
        ylabel(var_l2(vascao.factors{i}.ind))
    else
        loadings(vascao.factors{i}, 'Title','Case/control','VarsLabel', var_l2(sort(vascao.factors{i}.ind)));
    end

end


% Time Factor

i=2;
if isfield(vascao.factors{i},'scores')

    TvASCA = table((1:length(vascao.factors{i}.ind))',var_l2(sort(vascao.factors{i}.ind))', vascao.p(sort(vascao.factors{i}.ind),i),'VariableNames', {'Order','Label',sprintf('PvalueF%i',i)})
    
    vascao.factors{i}.loads = -vascao.factors{i}.loads;
    vascao.factors{i}.scoresV = -vascao.factors{i}.scoresV;

    scores(vascao.factors{i}, 'Title', 'Time', 'ObsClass',vascao.design(:,i));
    legend('T3', 'Labour');

    if length(vascao.factors{i}.ind)==1
        ylabel(var_l2(vascao.factors{i}.ind))
    else
        loadings(vascao.factors{i}, 'Title','Time', 'VarsLabel', var_l2(sort(vascao.factors{i}.ind))); 
    end
    
end


% Factor Individual
    
i=3;
if isfield(vascao.factors{i},'scores')

    TvASCA = table((1:length(vascao.factors{i}.ind))',var_l2(sort(vascao.factors{i}.ind))', vascao.p(sort(vascao.factors{i}.ind),i),'VariableNames', {'Order','Label',sprintf('PvalueF%i',i)})
    
    if vascao.factors{i}.stasig
        varPca(vascao.factors{i}.matrix);
        mspcPca(vascao.factors{i}.matrix, 'PCs', 1:2, 'ObsTest', vascao.factors{i}.matrix+vascao.residuals, 'PlotCal', false, 'ObsLabel', vascao.design(:,3), 'ObsClass', vascao.design(:,1), 'LimType', 1);
        title(sprintf('Factor %d',i));
        legend('Control', 'LFG');

    end
end


% 'Case/control' x 'Time'

i=1;
if isfield(vascao.interactions{i},'scores')

    TvASCA = table((1:length(vascao.interactions{i}.ind))',var_l2(sort(vascao.interactions{i}.ind))', vascao.p(sort(vascao.interactions{i}.ind),3+i),'VariableNames', {'Order','Label',sprintf('PvalueI%i',i)})
    
    scores(vascao.interactions{i}, 'Title', 'Int LFR/Ctrl x Time', 'ObsClass', vascao.design(:,1)*10+vascao.design(:,2));
    legend('Ctrl+T3','Case+T3','Ctrl+Labour','Case+Labour');

    if length(vascao.interactions{i}.ind)==1
        ylabel(var_l2(vascao.interactions{i}.ind))
    else
        loadings(model,'Title', 'Int LFR/Ctrl x Time', 'VarsLabel', var_l2(sort(vascao.interactions{i}.ind))); 
    end

end
