model = model_c;
sorting = true;


PlotVars.X = linspace(1,length(model.Coefficients.Estimate),...
    length(model.Coefficients.Estimate))';
PlotVars.Y = model.Coefficients.Estimate;
PlotVars.interval = coefCI(model);
% Sortieren der Faktoren nach absteigendem Betrag des Koeffizienten
if sorting
    [~, ind] = sort(abs(PlotVars.Y), 1, 'descend');
    PlotVars.Y = PlotVars.Y(ind);
    PlotVars.interval = PlotVars.interval(ind,:);
end
% Plot confidence intervals for the coefficients of the model
figure;
bar(PlotVars.Y)
hold on
errorbar(PlotVars.X, PlotVars.Y, PlotVars.interval(:,1), ...
    PlotVars.interval(:,2),'xr')
xlabel('Koeffizienten')
set(gca,'XTicklabel', model.CoefficientNames(ind))
hold off

clear PlotVars model ind sorting