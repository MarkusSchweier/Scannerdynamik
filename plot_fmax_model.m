%% Skript zum Plotten des Grenzfrequenzmodells gegen die Datenpunkte
clear all, close all, clc
% Laden des Modells
load([pwd, '\Daten\model_fmax.mat']);
% Extrahieren der Datenpunkte
[AEin, t, fmax] = dataset2vec(model_fmax.Variables);
AEin = cell2mat(AEin); t = cell2mat(t); fmax = cell2mat(fmax);
% Herausfiltern der Ausreißer
exc = dataset2vec(model_fmax.ObservationInfo(:,2));
exc = cell2mat(exc);
AEin = AEin(exc ~= 1); t = t(exc ~= 1); fmax = fmax(exc ~=1);
%Ableiten der konvexen Hülle
hullind = convhull(AEin, t);
% Berechenen der Modellwerte
    % Anlegen des Grids
    Amod = linspace(min(AEin), max(AEin), 100);
    tmod = linspace(min(t), max(t), 100);
    [Amod, tmod] = meshgrid(Amod, tmod);
    % Herausfiltern der Werte außerhalb des Versuchsraums
    out = ~inpolygon(Amod, tmod, AEin(hullind), t(hullind));
    Amod(out) = NaN; tmod(out) = NaN;
    % Vorhersage
    fmod = feval(model_fmax, Amod, tmod);
% Plotten
fig = figure;
set(fig,'Units','centimeters','Position',[2 5 7.75 7.75]);
    % Datenpunkte
    sg = scatter3(AEin, t, fmax);
    set(sg, 'Marker', '.', 'MarkerEdgeColor', 'k');
    sa = gca;
    set(sa, 'LineWidth', 1);
    xlab = xlabel(sa, 'Amplitude A_{Ein}', 'Units', 'normalized');
    ylab = ylabel(sa, 'Zeit t', 'Units', 'normalized');
    zlab = zlabel(sa, 'Grenzfrequenz f_{Grenz}', 'Units', 'normalized');
    set (xlab, 'FontSize', 13, 'FontName', 'Arial');
    set (ylab, 'FontSize', 13, 'FontName', 'Arial');
    set (zlab, 'FontSize', 13, 'FontName', 'Arial');
    set(sa, 'YDir', 'reverse', 'XDir', 'reverse');
    set(sa, 'xlim', [3.75 6], 'ylim', [0 1.75], 'zlim', [0 1200]);
    set(sa, 'XTick', [4:0.5:6], 'YTick', [0:0.5:1.5], 'ZTick', [0:300:1200]);
    XTL = {'4,0', '4,5', '5,0', 'mm', '5,5'};
    YTL = {'0,0', '0,5', 's', '1,5'};
    ZTL = {'0,0', '0,3', '0,6', 'kHz', '1,2'};
    set(sa, 'XTickLabel', XTL, 'YTickLabel', YTL, 'ZTickLabel', ZTL);
    set(sa, 'FontSize', 13, 'FontName', 'Arial');
    % Versuchsraumeinschränkung
    hold on;
    lp = plot3(sa, AEin(hullind), t(hullind), zeros(length(hullind)));
    set(lp, 'Color', 'k', 'LineWidth', 0.5);
    % Modell
    sp = surf(Amod, tmod, fmod);
    set(sp, 'LineStyle', 'none', 'FaceAlpha', 0.5);
        % Colormap
        cmap = [linspace(1, 0, 64); linspace(1, 0, 64); linspace(1, 0, 64)];
        cmap = cmap';
    colormap(cmap);
    hold off;