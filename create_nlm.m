function model = create_nlm(ds, modelfun, beta0)
% Funktion zur Erstellung nicht-linearer Modelle für alle im WS befindlichen
% Datasets
% ds:       Dataset-Array aus welchem das Modell erstellt werden soll
% modelfun: Function-Handle zur Modellfunktion @(b,x)
    % b:    nx1 Vektor von Koeffizienten
    % x:    Vektor der Faktoren
% beta0:    nx1 Vektor mit Startwerten für die Koeffizienten

% Erzeugung des NLM durch Regression
fprintf('-----\nNicht-Lineare-Regression des Dataset-Arrays.\n');
model = NonLinearModel.fit(ds, modelfun, beta0)

fprintf('Gütemaße\nR2:\t\t%.4f\nR2adj:\t%.4f\nAIC:\t%.4f\nBIC:\t%.4f\n', ...
        model.Rsquared.Ordinary,model.Rsquared.Adjusted,model.ModelCriterion.AIC,model.ModelCriterion.BIC);
end