function models = create_lm(mtype, p, l, u, Verbose)
% Funktion zur Erstellung linearer Modelle für alle im WS befindlichen
% Datasets
% p: Grenzwert für die Irrtumswahrscheinlichkeit
% l: Untergrenze des Models
% u: Obergrenze des Modells
% mtype: Modelltyp
% Verbose: Displayeinstellungen

% Auslesen der WS-Variablen
all_vars = evalin('base','whos');
% Identifikation aller Dataset-Arrays
models = struct;
for i = 1 : 1 : numel(all_vars)
    if isequal(all_vars(i).class, 'dataset')
        % Erzeugung des LM durch Regression
        fprintf('-----\nANOVA für das Dataset-Array: %s\n', all_vars(i).name);
        temp = LinearModel.stepwise(evalin('base', all_vars(i).name), mtype, 'Lower', l, 'Upper', u, 'Verbose', Verbose, 'Penter', p, 'PRemove', p + 1e-6);
        models = setfield(models, all_vars(i).name, temp);
        disp(temp.Coefficients);
        fprintf('Gütemaße\nR2:\t\t%.4f\nR2adj:\t%.4f\nAIC:\t%.4f\nBIC:\t%.4f\n', ...
                temp.Rsquared.Ordinary,temp.Rsquared.Adjusted,temp.ModelCriterion.AIC,temp.ModelCriterion.BIC);
    end
end
clear i j all_vars temp
end