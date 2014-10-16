function fmax = eval_fmax(A)
% Funktion zur Berechnung der Grenzfrequenz für eine Dauerhafte Oszillation
% A: Soll-Amplitude auf dem Bauteil in mm
% fmax: Maximal-Frequenz in Hz

% Laden der Daten
load([pwd, '\Daten\data_fmax.mat']);
if A < min(Avec)
    fmax = fvec(1);
elseif A >= min(Avec) && A <= max(Avec)
    fmax = interp1(Avec,fvec,A,'linear','extrap');
else
    fmax = fvec(length(fvec));
end
