function ds_coded = code_ds(ds, lb, ub)
% Funktion zum codieren eines Versuchsplans in einem Dataset-Array
% ds: Dataset-Array
% ds_coded: Codiertes Data-Set-Array
% lb: Untere Grenze der Kodierung
% ub: Obere Grenze der Kodierung

% Umwandeln des DS in ein Cell-Array
c = dataset2cell(ds);
% Entfernen der Kopfzeile
c = c(2:end,:);
% Umwandeln in eine Matrix
m = cell2mat(c);
% Aufteilen in Versuchsplan-Matrix und Ergebnis-Vektor
% Annahme: Nur letzte Spalte enthält Ergebnisse
X = m(:,1:end-1);
Y = m(:,end);
% Kodieren der Faktoreinstellungen
for i = 1:1:length(X(1,:))
    % Auslesen Maximal-Wert
    maxim = max(X(:,i));
    % Auslesen Minimal-Wert
    minim = min(X(:,i));
    % Kodieren
    Xc(:,i) = 1 / (maxim - minim) * ((ub - lb) * X(:,i) + lb * maxim - ub * minim);
end
% Zusammensetzten des neuen DS
ds_coded = replacedata(ds,horzcat(Xc,Y));
   
end