%% Ferrer et al. INTESTINAL MICROBIOTA INFLAMMATORY PROFILE IN LATE-FETAL GROWTH RESTRICTION WOMEN, 2025. 
% Script to read and prepare the code.
%
% Software preparation: Install MEDA-Toolbox v1.8
%
% coded by: Jose Camacho (josecamacho@ugr.es)
%       Óliver Aleksandrei Polushkina (oap327@gmail.com)
% last modification: 7/Abr/2025
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

%% Load the data

clear
close all
clc

data = importdata('BD_ESTADISTICA_PI1701215_CGLL.xlsx');

X = data.data.BD_analisis(1:95, 3:end);

obs_l = data.textdata.BD_analisis(3:end-1, 2);
var_l = data.textdata.BD_analisis(2, 3:end);

class(find(ismember(obs_l, 'control'))) = 1;
class(find(ismember(obs_l, 'CIR'))) = 2;
class(find(ismember(obs_l, 'PEG'))) = 3;
class = class';
class_l = data.textdata.BD_analisis(2, 2);

save data_base.mat X var_l obs_l class


%% Blood-hemodynamic data preparation

clear
close all
clc
load data_base.mat

X_blood = X(1:95, 1:12);
X_hemo = X(1:95, 14:25);

% Etiquetas variables
var_blood = var_l(1:end, 1:12);
var_hemo = var_l(1:end, 14:25);

var_final = {'Hemoglobin', 'Hematocrite', 'MCV', 'Platelets', 'Systolic BP', ...
    'Diastolic BP', 'Mother weight', 'BMI'};

% Joining matrices and removing PEG individuals
X_bh = [X_blood, X_hemo];
var_bh = [var_blood, var_hemo];
ind = find(class == 3);
X_bh(ind, :) = [];
class(ind) = [];

% Rearrenge matrix
ind_T1 = [find(contains(var_bh, '1T')), find(contains(var_bh, '12'))];
ind_T2 = [find(contains(var_bh, '2T')), find(contains(var_bh, '24'))];
ind_T3 = [find(contains(var_bh, '3T')), find(contains(var_bh, '32'))];

X_tiempo = [];
for i = 1:length(ind_T1)
    X_tiempo(:, i) = [X_bh(:, ind_T1(i)); X_bh(:, ind_T2(i)); X_bh(:, ind_T3(i))];
end

tiempo = repelem(1:3, length(class))';
id = repmat(1:length(class), 1, 3);
F_ASCA = [[class; class; class], tiempo, id'];
X_ASCAt = X_tiempo;

save blood_hemo_vasca.mat X_ASCAt F_ASCA id var_bh var_final


%% Inflammatory Biomark dataprep

clear
close all
clc
load data_base.mat

% Data matrix + labels
X_biomark = X(1:95, 34:75);
var_biomark = var_l(1:end, 34:75);

% Remove PEG individuals
ind = find(class == 3);
X_biomark(ind, :) = [];
class(ind, :) = [];

var_final = {'RARRESS', 'LBP', 'IFN-G', 'IL-10', 'IL-23','IL-6', ...
    'IL-8', 'Adiponectin', 'Resistin', 'IL-15', 'VEGF', 'Leptin', ...
    'Adipo/Lep', 'Adipo/Resis'};

id = repmat(1:length(class), 1, 3);
F_mom_baby = [[class; class; class], repelem(1:3, length(class))', id'];
ind_pre = find(contains(var_biomark, '_pre'));
ind_parto = find(contains(var_biomark, '_parto'));
ind_stage = find(contains(var_biomark, 'fetal'));

X_tiempo = [];
for i = 1:14
    X_tiempo(:, i) = [X_biomark(:, ind_pre(i)); X_biomark(:, ind_parto(i)); X_biomark(:, ind_stage(i))];
end

X_ASCAt = X_tiempo;
F_ASCA = F_mom_baby;

save biomark_vasca.mat X_ASCAt F_ASCA id var_biomark var_final


%% Bacteria dataprep

clear
close all
clc
load data_base.mat

% Data matrix of bacteria and labels
X_bact = X(:, 76:end);
var_bact = var_l(:, 76:end);

% Removing PEG individuals
ind = find(class == 3);
X_bact(ind, :) = [];
class(ind) = [];

% Create F
id = repmat(1:length(class), 1, 4);
F = [repmat(class, 4 ,1), repelem(1:4, length(class))', id'];

% Rearrange data
var_nombre = regexprep(var_bact, '3R', '1');
var_nombre = regexprep(var_nombre, 'MADRE', '2');
var_nombre = regexprep(var_nombre, 'MECONIO', '3');
var_nombre = regexprep(var_nombre, 'CUARENTENA', '4');

[var_ord, ind] = sort(var_nombre);
X_tiempo = X_bact(:, ind);
var_final = unique(regexprep(var_ord, '[\d"]', ''));

for i=1:length(var_final), var_final{i} = strrep(var_final{i},'__','-'); end
for i=1:length(var_final), var_final{i} = strrep(var_final{i},'_',''); end

tmp1 = find(contains(var_ord, '_1'));
tmp2 = find(contains(var_ord, '_2'));
tmp3 = find(contains(var_ord, '_3'));
tmp4 = find(contains(var_ord, '_4'));

% Data matrix, each column a bacteria
X_final = [X_tiempo(:, tmp1); X_tiempo(:, tmp2); X_tiempo(:, tmp3); X_tiempo(:, tmp4)];

X_ASCAt = X_final;
F_ASCA = F;

save bacteria_vasca.mat X_ASCAt F_ASCA var_final var_bact