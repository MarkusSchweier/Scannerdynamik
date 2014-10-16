function [varargout] = dataset2vec(ds)
% Funktion zur Konvertierung eines Dataset-Arrays in einzelne Vektoren des Datentyps Cell
% ds:       Dataset-Array mit N-1 Faktoren und einer Zielgröße an letzter Stelle
% varagout: N-1 Faktoren und Zielgröße als einzelne Vektoren

% Konvertieren des ds in eine Cell-Array
ca = table2cell(ds);
% Auslesen der Anzahl an Rückgabevektoren
num = numel(ca(1,:));
% Prüfen, ob die richtige Zahl an Ausgabevariablen angefordert wurde
if nargout ~= num
    error('Falsche Anzahl an Rückgabewerten beim Fuktionsaufruf gefordert.');
end
% Entfernen der Namen aus ca
ca = ca(2:end,:);
% Auslesen der Vektoren für die Faktoren
for i = 1 : 1 : num
    varargout{i} = ca(:,i);
end
end