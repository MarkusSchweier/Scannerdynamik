function [varargout] = dataset2vec(ds)
% Funktion zur Konvertierung eines Dataset-Arrays in einzelne Vektoren des Datentyps Cell
% ds:       Dataset-Array mit N-1 Faktoren und einer Zielgr��e an letzter Stelle
% varagout: N-1 Faktoren und Zielgr��e als einzelne Vektoren

% Konvertieren des ds in eine Cell-Array
ca = table2cell(ds);
% Auslesen der Anzahl an R�ckgabevektoren
num = numel(ca(1,:));
% Pr�fen, ob die richtige Zahl an Ausgabevariablen angefordert wurde
if nargout ~= num
    error('Falsche Anzahl an R�ckgabewerten beim Fuktionsaufruf gefordert.');
end
% Entfernen der Namen aus ca
ca = ca(2:end,:);
% Auslesen der Vektoren f�r die Faktoren
for i = 1 : 1 : num
    varargout{i} = ca(:,i);
end
end