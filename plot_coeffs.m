% Festlegung des Models und der Koeffizientenreihenfolge
mod = 'A_Ist_c';
sorting = true;
% Größeneinstellungen für den Plot
B = 7.75; % Breite in cm
H = 5; % Höhe in cm
% Hinzufügen des Datenpfads
addpath([pwd,'/Daten']);
load('lm.mat');
% Auslesen der Daten
model = eval(['models.', mod]);
PlotVars.X = linspace(1,length(model.Coefficients.Estimate),...
    length(model.Coefficients.Estimate))';
PlotVars.Y = model.Coefficients.Estimate;
% Schätzen der Konfidenzintervalle
PlotVars.interval = coefCI(model);
% Sortieren der Faktoren nach absteigendem Betrag des Koeffizienten
if sorting
    [~, ind] = sort(abs(PlotVars.Y), 1, 'descend');
    PlotVars.Y = PlotVars.Y(ind);
    PlotVars.interval = PlotVars.interval(ind,:);
end
% Plotten der Koeffizienten
figure;
set(gcf,'Units','centimeters');
set(gcf,'ActivePositionProperty','Position');
set(gcf,'Position',[5 5 B H]);
bar(PlotVars.Y,'FaceColor',[1 1 1] .* 221 / 255)
% Grid
set(gca,'XGrid','off','YGrid','on','ZGrid','off')
% Achslimit
xlim([0, length(PlotVars.X) + 1]);
hold on
% Berechnen der CI-Grenzen
u = PlotVars.interval(:,1) - PlotVars.Y;
l = PlotVars.interval(:,2) - PlotVars.Y;
% Plotten der CIs
errorbar(PlotVars.X, PlotVars.Y, u, l,'.k');
% Achsbeschriftungen
set(gca, 'LineWidth', 1.0);
set(gca, 'FontSize', 13);
set(gca, 'Fontname', 'Arial');
xlabel('Faktoren')
ylabel('Koeffizient \beta');
set(gca,'XTicklabel', model.CoefficientNames(ind))
% Größe
set(gca, 'ActivePositionProperty', 'Outerposition');
set(gca, 'Units', 'normalized');
set(gca, 'Outerposition', [0 0 1 1]);
hold off
% Formatieren der XTickLabel
load(['ticks_', mod]);
[~,~] = format_ticks(gca,ticks,[],[],[],[],[],[]);
% Löschen der temporären Variablen
clear PlotVars model ind sorting u l mod h ticks B H