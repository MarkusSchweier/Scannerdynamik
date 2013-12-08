% Skript zur Erzeugung des Ergebnisplots für die Grenzfrequenzen
clear all, close all, clc
% Laden und Vorbereiten der Daten
    % Laden
    load([pwd,'\Daten\model_fmax.mat']);
    [Adata, tdata, fdata] = dataset2vec(model_fmax.Variables);
    excdata = dataset2vec(model_fmax.ObservationInfo(:,2));
    Adata = cell2mat(Adata); tdata = cell2mat(tdata); fdata = cell2mat(fdata); excdata = cell2mat(excdata);
    % Entfernen der Ausreisen
    Adata = Adata(excdata ~= 1); tdata = tdata(excdata ~= 1); fdata = fdata(excdata ~= 1);
    % Erzeugen der konvexen Hülle
    hullind = convhull(Adata, tdata);
% Berechnung der Modellvorhersagen
    % Anlegen eines Grids
    Amodel = linspace(min(Adata), max(Adata), 100);
    tmodel = linspace(min(tdata), max(tdata), 100);
    [x1, x2] = meshgrid(Amodel, tmodel);
    IN = inpolygon(x1, x2, Adata(hullind), tdata(hullind));
    x1(~IN) = NaN; x2(~IN) = NaN;
    fmodel = feval(model_fmax, x1, x2);
    
% Plotten
fig = figure;
    % Daten
    scatter3(Adata, tdata, fdata, 'k', '.');
    set(gca, 'XDir', 'reverse', 'YDir', 'reverse', 'LineWidth', 1);
    % Hülle
    hold on;
    plot3(Adata(hullind), tdata(hullind), zeros(length(hullind)), 'LineWidth', 0.5, 'Color', [0 0 0]);
    % Modelldaten
    surf(x1, x2, fmodel, 'LineStyle', 'none');
    alpha(0.7);
