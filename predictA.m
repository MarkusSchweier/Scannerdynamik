function [AEin, fmax, fvec, boundsAEin, boundsfmax] = predictA(f, ASoll, t)
% Funktion zur Vorhersage der Eistellamplitude unter Ber�cksichtigung des
% D�mpfungsmodells und des Modells der Leistungsaufnahmegrenze
% AEin:         n x p Matrix der Einstellamplituden
% fmax:         n x p Matrix mit Grenzfrequenzen f�r Schwei�zeit t
% fvec:         1 x p Vektor mit den zur Schwei�zeit zugeh�rigen Grenzfrequenzen f�r die gegebenen ASoll
% boundsAEin:   n x p Matrix mit logischen Werten, die angibt, ob die AEin-Prognose mit Werten gesch�tzt wurde, die im Definitionsbereich des Modells liegen
% boundsfmax:   n x p Matrix mit logischen Werten, die angibt, ob die fmax-Prognose mit Werten gesch�tzt wurde, die im Definitionsbereich des Modells liegen
% f:            n x 1 Vektor mit den gew�nschten Frequenzen
% ASoll:        1 x p Vektor mit den gew�nschten Amplituden auf dem Bauteil
% t:            Skalar der Schwei�zeit, f�r welche die Leistungsaufnahmegrenze gelten soll

%% Pr�fen und ggf. Anpassen der �bergebenen Variablen
% Anzahl der Argumente
    if nargin ~= 3
        error('Falsche Anzahl von Eingabewerte.');
    end
% Datentyp der Argumente
    if ~isvector(f) || ~isvector(ASoll) || ~isscalar(t)
        error(['Falscher Datentyp mindestens eines Argumentes.\n', ...
        'f und ASoll m�ssen Vektoren sein. t muss ein Skalar sein.']);
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
%% Erzeugen von Spaltenvektoren des Bereichs f�r die vektorisierte Verarbeitung
[x2, x1] = meshgrid(ASoll, f);
%% Festlegung der Modelldaten
% D�mpfungsmodell
b = 1129.5;
Arange = [0.1 5.5]; %[0.105 5.530];
frange = [100, 2000]; %[101, 1999];
% Leistungsaufnahmegrenze-Modell
    % Laden
    load([pwd,'\Daten\model_fmax.mat']);
    % Erzeugen der konvexen H�lle
    [Adata, tdata, ~] = dataset2vec(model_fmax.Variables);
    excdata = dataset2vec(model_fmax.ObservationInfo(:,2));
    Adata = cell2mat(Adata); tdata = cell2mat(tdata); excdata = cell2mat(excdata);
    Adata = Adata(excdata ~= 1); tdata = tdata(excdata ~= 1);
    hullind = convhull(Adata, tdata);
    trange = model_fmax.VariableInfo.Range{2,1};
%% Berechnung der R�ckgabewerte
AEin = x2 .* exp((x1./b).^2);
fmax = feval(model_fmax, AEin, repmat(t,size(AEin)));
% Suchen der Elemente au�erhalb des Definitionsbereichs
    boundsAEin = false(size(AEin));
    boundsfmax = false(size(AEin));
    boundsAEin(x1 < min(frange) | x1 > max(frange) | x2 < min(Arange) | x2 > max(Arange)) = 1;
    boundsfmax = ~inpolygon(AEin,repmat(t,size(AEin)),Adata(hullind),tdata(hullind));
    boundsfmax = boundsfmax | boundsAEin;
% Reduzierung der Ergebnismatrizen um nicht valide prognostizierte Daten
if sum(sum(boundsAEin)) > 0
    fprintf(['-----\nEs kann nicht f�r alle Faktorkombinationen eine Einstellamplitude vorhergesagt werden.\n', ...
    'Die entsprechenden Elemente werden zu NaN gesetzt.\n-----\n']);
    AEin(boundsAEin) = NaN;
end
if sum(sum(boundsfmax)) > 0
    fprintf(['-----\nEs kann nicht f�r alle Faktorkombinationen eine Grenzfrequenz vorhergesagt werden.\n', ...
    'Die entsprechenden Elemente werden zu NaN gesetzt.\n-----\n']);
    fmax(boundsfmax) = NaN;
end
% �berpr�fung der Leistungsaufnahmegrenze und Auslesen eines f_Grenz-Vektors der Gr��e n x 1
temp = false(size(x1));
temp(bsxfun(@lt,fmax,x1))=1;
fvec = NaN(1,length(ASoll));
% Schleife durch alle Spalten
for i = 1 : 1 : length(temp(1,:))
    % Finden der untersten 1
    row = find(temp(:,i) == 1);
    % Setzten aller Werte dar�ber auf 1
    temp(min(row):end,i) = 1;
    % Auslesen der zugeh�rigen Grenzfrequenz
    if ~isempty(row)
        fvec(i) = f(min(row));
    end
end
AEin(temp) = NaN;
end