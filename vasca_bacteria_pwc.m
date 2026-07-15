%% Ferrer et al. INTESTINAL MICROBIOTA INFLAMMATORY PROFILE IN LATE-FETAL GROWTH RESTRICTION WOMEN, 2025. 
% Script to analyze bacterial data
%
% We consider the following model: 
%
%   X = A + B + C(A) + AB + E
%
% with A: Class (Control vs IUGR), B: Time (T3, Labour), C(A): Individual
%
% Software preparation: Install MEDA-Toolbox v1.14
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 13/Jul/2026
%
% Copyright (C) 2026  University of Granada, Granada
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

load bacteria_vasca.mat

X = X_ASCAt;
Xraw = X;
F = F_ASCA;
for i=1:length(var_final), var_l{i} = char(var_final(i)); end

ind = find(F(:,2)==3);% do not consider newborns
X(ind,:) = [];
Xraw(ind,:) = [];
F(ind,:) = [];

ind = find(F(:,2)==4);% do not consider newborns
X(ind,:) = [];
Xraw(ind,:) = [];
F(ind,:) = [];


clear X_ASCAt F_ASCA  

%% PWC v1.8

addpath(genpath('MEDA-Toolbox-1.8\toolbox'))

which('powercurveVSold.m')
which('parglmVS.m')
which('parglmMC.m')

D.N = size(F,1);
D.M = size(X,2);
D.k = [0.1 0.3 1 0.1];        

[PCmeanA, PCrepA, powercurveoA] = powercurveVSold(D, F, 'Alpha', pvalue, 'Repetitions', 200, 'Permutations', 200, 'Model', [1 2], 'Nested', [1 3], 'Random', [0 0 1]);
ax = gca;
box(ax, 'on');
% ax.XLabel.FontSize = 32;
% ax.YLabel.FontSize = 32;
% ax.FontSize = 20;
ylim(ax, [0 1]);
title('');
legend({'F1: LFR/Ctrl','F2: Time','F3: Individual','Int Class x Time'})


save bacteria_vasca_pwc

%% PWC v1.14

addpath(genpath('D:\Curro\Software\Datharsis\DathaTest\toolbox'))

which('powercurveVS.m')
which('parglmVS.m')
which('parglmMC.m')

D.N = size(F,1);
D.M = size(X,2);
D.k = [0.1 0.3 1 0.1];        

[PCmeanA14, PCrepA14, powercurveoA14] = powercurveVS(D, F, 'Alpha', pvalue, 'Repetitions', 200, 'Permutations', 200, 'Model', [1 2], 'Nested', [1 3], 'Random', [0 0 1]);
ax = gca;
box(ax, 'on');
% ax.XLabel.FontSize = 32;
% ax.YLabel.FontSize = 32;
% ax.FontSize = 20;
ylim(ax, [0 1]);
title('');
legend({'F1: LFR/Ctrl','F2: Time','F3: Individual','Int Class x Time'})


save bacteria_vasca_pwc PCmeanA14 PCrepA14 powercurveoA14 -APPEND