function [AEin, fmax, boundsAEin, boundsfmax] = predictA(f, ASoll, t)
% Funktion zur Vorhersage der Eistellamplitude unter Berücksichtigung des
% Dämpfungsmodells und des Modells der Leistungsaufnahmegrenze
% AEin:         n x p Matrix der Einstellamplituden
% fmax:         n x p Matrix mit Grenzfrequenzen für Schweißzeit t
% boundsAEin:   n x p Matrix mit logischen Werten, die angibt, ob die AEin-Prognose mit Werten geschätzt wurde, die im Definitionsbereich des Modells liegen
% boundsfmax:   n x p Matrix mit logischen Werten, die angibt, ob die fmax-Prognose mit Werten geschätzt wurde, die im Definitionsbereich des Modells liegen
% f:            n x 1 Vektor mit den gewünschten Frequenzen
% ASoll:        1 x p Vektor mit den gewünschten Amplituden auf dem Bauteil
% t:            Skalar der Schweißzeit, für welche die Leistungsaufnahmegrenze gelten soll

%% Prüfen und ggf. Anpassen der übergebenen Variablen
% Anzahl der Argumente
    if nargin ~= 3
        error('Falsche Anzahl von Eingabewerte.');
    end
% Datentyp der Argumente
    if ~isvector(f) || ~isvector(ASoll) || ~isscalar(t)
        error(['Falscher Datentyp mindestens eines Argumentes.\n', ...
        'f und ASoll müssen Vektoren sein. t muss ein Skalar sein.']);
    end
% Ggf. transponieren
    sizef = size(f);
    sizeASoll = size(ASoll);
    if sizef(2) ~= 1
        f = f';
    end
    if sizeASoll(1) ~= 1
        ASoll = ASoll';
    end
%% Erzeugen von Spaltenvektoren des Bereichs für die vektorisierte Verarbeitung
[x2, x1] = meshgrid(ASoll, f);
%% Festlegung der Modelldaten
% Dämpfungsmodell
b = 1129.5;
Arange = [0.105 5.530];
frange = [101, 1999];
% Leistungsaufnahmegrenze-Modell
    % Laden
    load([pwd,'\Daten\model_fmax.mat']);
    % Erzeugen der konvexen Hülle
    [Adata, tdata, ~] = dataset2vec(model_fmax.Variables);
    excdata = dataset2vec(model_fmax.ObservationInfo(:,2));
    Adata = cell2mat(Adata); tdata = cell2mat(tdata); excdata = cell2mat(excdata);
    Adata = Adata(excdata ~= 1); tdata = tdata(excdata ~= 1);
    hullind = convhull(Adata, tdata);
    trange = [0.1412 1.5];
%% Berechnung
AEin = x2 .* exp((x1./b).^2);
fmax = feval(model_fmax, AEin, repmat(t,size(AEin))); % fmax = b1 + b2./AEin + b3/t + b4./AEin.^2 + b5/t^2 + b6./(AEin*t);
% Suchen der Elemente außerhalb des Definitionsbereichs
    boundsAEin = false(size(AEin));
    boundsfmax = false(size(AEin));
    boundsAEin(x1 < min(frange) | x1 > max(frange) | x2 < min(Arange) | x2 > max(Arange)) = 1;
    % Aufstellen boundsAEin
%     [r1, c1, ~] = find(x1 < min(frange) | x1 > max(frange));
%     [r2, c2, ~] = find(x2 < min(Arange) | x2 > max(Arange));
%     r = [r1; r2]; c = [c1; c2];
%     boundsAEin(sub2ind(size(boundsAEin), r', c')) = 1;
if t < min(trange) || t > max(trange)
    boundsfmax(:,:) = 1;
else
    boundsfmax = ~inpolygon(AEin,repmat(t,size(AEin)),Adata(hullind),tdata(hullind));
    boundsfmax = boundsfmax | boundsAEin;
end
% Reduzierung der Ergebnismatrizen um nicht valide prognostizierte Daten
if sum(sum(boundsAEin)) > 0
    fprintf(['-----\nEs kann nicht für alle Faktorkombinationen eine Einstellamplitude vorhergesagt werden.\n', ...
    'Die entsprechenden Elemente werden zu NaN gesetzt.\n-----\n']);
    AEin(boundsAEin) = NaN;
end
if sum(sum(boundsfmax)) > 0
    fprintf(['-----\nEs kann nicht für alle Faktorkombinationen eine Grenzfrequenz vorhergesagt werden.\n', ...
    'Die entsprechenden Elemente werden zu NaN gesetzt.\n-----\n']);
    fmax(boundsfmax) = NaN;
end
% Überprüfung der Leistungsaufnahmegrenze
temp = false(size(x1));
temp(bsxfun(@lt,fmax,x1))=1;
% Schleife durch alle Spalten
for i = 1 : 1 : length(temp(1,:))
    % Finden der untersten 1
    row = find(temp(:,i) == 1);
    % Setzten aller Werte darüber auf 1
    temp(min(row):end,i) = 1;
end
AEin(temp) = NaN;
end