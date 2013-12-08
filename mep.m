function [] = mep(mod, B, H, alpha)
% Funktion zur Ausgabe eines MainEffectsPlot
% mod:  Regressionsmodell der Klasse LinearModel oder NonLinearModel
% B:    Breite der Figure
% H:    Höhe der Figure
% alpha:  1-alpha CIs werden geplottet; falls alpha = 0, werden keine CIs geplottet

% Einstellwerte
num = 100; % Anzahl der diskreten Funktionswerte für einen MEP
part = 0.05; % Prozentualer Abstand zwischen zwei Subplots

% Auslesen der im Modell enthaltenen Haupteffekte und deren Range
effects.names = cell(1);
effects.ranges = [0 0];
z = 1;
for i = 1 : 1 : length(mod.VariableInfo) - 1
    if mod.VariableInfo.InModel(i,1)
        effects.names{z,1} = mod.VariableNames{i,1};
        effects.ranges(z,:) = mod.VariableInfo{i,2};
        effects.num(z,1) = i; 
        z = z + 1;
    end
end
% Berechnung der Zielgrößenwerte und CIs
% Anlegen der Vektoren
x = zeros(num, mod.NumVariables - 1, numel(effects.names));
y = zeros(num,1,numel(effects.names));
ci = zeros(num,2,numel(effects.names));
for i = 1 : 1 : numel(effects.names)
    for j = 1 : 1 : mod.NumVariables - 1 
        if mod.VariableInfo{j, 3} && isequal(effects.names{i,1}, mod.VariableNames{j,1})
            x(:,j,i) = linspace(effects.ranges(i,1), effects.ranges(i,2));
        elseif mod.VariableInfo{j, 3} && ~isequal(effects.names{i,1}, mod.VariableNames{j,1})
            x(:,j,i) = mean(effects.ranges(i,:));
        end
    end
    % Berechnung der y-Werte
    for k = 1 : 1 : num
        [y(k,1,i), ci(k,:,i)] = predict(mod,x(k,:,i));
    end
end

% Erzeugen der Figure
fig = figure;
set(fig, 'Units', 'centimeters', 'ActivePositionProperty', 'Position');
set (fig, 'Position', [5 5 B H]);
% Hinzufügen der Subplots
for i = 1 : 1 : numel(effects.names)
    sph(i) = subplot(1, numel(effects.names), i);
    set(sph(i), 'LineWidth', 1.0);
    set(sph(i), 'FontSize', 13);
    box on;
    grid on;
    hold on;
    plot(x(:,effects.num(i),i), y(:,1,i), 'LineWidth', 1.0, 'LineStyle', '-', 'Color', [0 0 0]);
    plot(x(:,effects.num(i),i), ci(:,1,i), 'LineWidth', 0.5, 'LineStyle', '-', 'Color', [0 0 0] .* 221 / 255);
    plot(x(:,effects.num(i),i), ci(:,2,i), 'LineWidth', 0.5, 'LineStyle', '-', 'Color', [0 0 0] .* 221 / 255);
    xlabel(effects.names{i,1});
    if i == 1
        ylabel(mod.ResponseName);
        ylims = ylim;
    else
        % Auslesen der ylims
        ylims(1,1) = min(min(ylims), min(get(sph(i), 'ylim')));
        ylims(1,2) =  max(max(ylims), max(get(sph(i), 'ylim')));
        % Entfernen der Ticks
        set(sph(i), 'yTickLabel',[]);
    end
    hold off;  
end
% Anpassen der Größen der Ades-Objekte und der ylim
% Verschieben des letzten Subplots an den rechten Rand
xfirst = get(sph(1),'OuterPosition');
xfirst(1) = 0;
set(sph(1), 'Outerposition', xfirst);
xlast = get(sph(length(sph)),'OuterPosition');
xlast(1) = 1 - xlast(3);
set(sph(length(sph)), 'Outerposition', xlast);
xmin = get(sph(1),'Position');
xmin = xmin(1);
xmax = get(sph(length(sph)),'Position');
xmax = xmax(1) + xmax(3);
b = ((xmax - xmin) - (length(sph) - 1) * part) / length(sph);
h = get(sph(1), 'Position');
for i = 1 : 1 : numel(effects.names)
    % Übertragen der ylims
    set(sph(i), 'YLimMode', 'manual', 'YLim', ylims);
    % Setzen der Position
    set(sph(i), 'ActivePositionProperty', 'position');
    pos = get(sph(i),'Position');
    pos(3) = b;
    x1 = get(sph(1), 'position');
    x1 = x1(1);
    pos(1) = x1 + (i-1) * (b + part);
    set(sph(i), 'position', pos);
    if h(2) < pos(2)
        h(2) = pos(2);
    end
    if h(4) > pos(4)
        h(4) = pos(4);
    end
end
for i = 1 : 1 : numel(effects.names)
    % Einstellen der y-Werte
    pos = get(sph(i), 'Position');
    pos(2) = h(2);
    pos(4) = h(4);
    set(sph(i), 'Position', pos);
end
clear i j k z
end