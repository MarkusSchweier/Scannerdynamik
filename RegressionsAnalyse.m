% Plot der verschiedenen Zielgrößen und des VP in 2D/3D Scatter Plots
%% Wahl der Variablen und Modelle
clear all
close all
clc

% Abfrage: gefiltert oder ungefiltert
filter = input(['Welche Daten sollen verwendet werden? f(gefiltert)' ...
    'oder of(ohne Filter): '],'s');
clc

% Laden der Zielgrößen
fprintf('Welche der folgenden Zielgrößen soll geladen werden?\n\n')
fprintf(['  DirDis1_B                DirDis1_E\n  DirDis2_B         '...
    '       DirDis2_E\n'])
fprintf('  spatter_per_lengthB      spatter_per_lengthE\n')
fprintf('  VeloDis_B                VeloDis_E\n\n')
% Benutzerabfrage nach Zielgröße
ziel = input('Eingabe der Zielgröße: ', 's');
% Laden der gewünschten Daten: gefiltert oder ohne Filter
if (strcmp('f', filter))
    response = load(['..\4.) Berechnen und Abspeichern der Zielgrößen\'...
        'Speicherung der Zielgrößen\gefiltert\' ziel '.mat']);
elseif (strcmp('of', filter))
    response = load(['..\4.) Berechnen und Abspeichern der Zielgrößen\'...
        'Speicherung der Zielgrößen\nicht gefiltert\' ziel '.mat']);
else
    fprintf('Bitte f oder of wählen.')
end
response = struct2cell(response);
response = cell2mat(response);

clc

% Laden des VP
VP = load('.\VP_uncoded.mat');
VP = VP.VP_uncoded;

load('.\group.mat')

% Wahl der VP-Größe
sizeVP = input('Sollen 50(5) oder 70(7) Punkte ausgewählt werden?: ', 's');

% Größe des Konfidenzintervalls
confi = 0.1;

% Skallierung der y-Achse für IntPlot
PlotVars.IntyScale = 5;

% Schrittweite für Skallierung der y-Achse für IntPlot
PlotVars.Step = 1;

% Schriftgröße
FontSize = 8;

% Initialisieren des Cell für die Modelle
Models = cell(3,1);

% Wert für Lambda für Box-Cox-Transformation
lambda = 0.5;

% Wahl des auszuwertenden Modells
fprintf('\nmdl1 --> LinearesModel mit ANOVA\n')
fprintf('mdl2 --> LinearesModel mit vorgegebener T-Matrix + RobustOpts\n')
fprintf('mdl3 --> NichtLinearesModel\n\n')
ModelChoice = input(['Welches Modell soll gewählt werden? mdl1(1),' ...
    ' mdl2(2) oder mdl3(3): '], 's');
ModelChoice = str2double(ModelChoice);

% Wahl ob Box-Cox-Transformiert werden soll oder nicht
Box_Cox = input(['\n\nSoll eine Box-Cox-Transformation mit Lambda='...
    num2str(lambda) ' durchgeführt werden?: Ja(1) oder Nein(2) '], 's');
Box_Cox = str2double(Box_Cox);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beschränkung auf 50 Versuche %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(sizeVP, '5'))
    VP = VP(1:end-20,:);
    response = response(1:end-20);
    group = group(1:end-20);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot
gplotmatrix(VP,response,group)

% Variable für Versuchsplan dataset
ds = cell(2,1);

% Box-Cox-Transformation
yTrans = boxcox(lambda, response);
ds{1,1} = mat2dataset([VP, yTrans], 'VarNames', ...
    {'a', 'f', 'v', 'p', ziel});

% Box-Cox-Back-Transformation
box_back = @(x) (lambda*(1/lambda + x)).^(1/lambda);

% Erstellen eines Datasets mit VP und Zielgröße
ds{2,1} = mat2dataset([VP, response], 'VarNames', ...
    {'a', 'f', 'v', 'p', ziel});
clc

% Erstellen eines linearen Regressionsmodells
Models{1,1} = LinearModel.stepwise(ds{Box_Cox,1}, 'quadratic', 'Lower', ...
    'constant', 'Upper', 'quadratic', 'Verbose', 1);

% Erstellen eines linearen Regressionsmodells
load('.\T.mat')
Models{2,1} = LinearModel.fit(ds{Box_Cox,1}, T, 'RobustOpts', 'fair', ...
    'ResponseVar', 5);

% create non-linear models
modelfun = @(b,x) b(1) + b(2)*x(:,1) + b(3)*x(:,2) + b(4)*x(:,3) + ...
    b(5)*x(:,4) + b(6)*(x(:,1).*x(:,4)) + b(7)*(x(:,1).^2);
beta0 = [0;0;0;0;0;0;0];
Models{3,1} = NonLinearModel.fit(ds{Box_Cox,1}, modelfun, beta0);

% Print R^2 value of used model
fprintf('\n\nR^2adj: %1.3f\n\n', Models{ModelChoice,1}.Rsquared.Adjusted)

%% Outlier Removal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Outlier entfernen mittels RESIDUEN
outl_single = [];
outl_top = [];
outl_bottom = [];
bad = input('Sollen outlier mit RESIDUEN entfernt werden? Ja(j) oder Nein(n) ', 's');
while (strcmp(bad, 'j'))
    % Identifizieren der outlier mit Cooks Distance and normal probability plot
    figure(1);
    plotDiagnostics(Models{ModelChoice,1},'cookd') %Cooks Distance
    grid on
    figure(2);
    plotResiduals(Models{ModelChoice,1},'probability')%probability plot
    grid on
    
    
    % Benutzerabfrage für outlier Identifizierung
    single = input('Einzel Versuch entfernen: ', 's');
    top = input('Residual für outlier cut-off TOP: ', 's');
    bottom = input('Residual für outlier cut-off BOTTOM: ', 's');
    
    % Umwandlung in double
    single = str2double(single);
    top = str2double(top);
    bottom = str2double(bottom);
    
    % Finden der outlier in Datensatz
    top = find(Models{ModelChoice,1}.Residuals.Raw > top);
    bottom = find(Models{ModelChoice,1}.Residuals.Raw < bottom);
    
    % anfügen an alte outlier
    if (~isnan(single))
        outl_single = cat(1, outl_single, single);
    end
    outl_top = cat(1, outl_top, top);
    outl_bottom = cat(1, outl_bottom, bottom);
    
    % Zusammenführen der angegebenen outlier
    outl = cat(1, outl_single, outl_top, outl_bottom);
    
    % Erstellen des neuen Modells ohne die outlier
    Models{1,1} = LinearModel.stepwise(ds{Box_Cox,1}, 'quadratic', ...
        'Upper', 'quadratic', 'Exclude', outl, 'Verbose', 0);
    
    % Erstellen des neuen linearen Modells mit T-Matrix
    Models{2,1} = LinearModel.fit(ds{Box_Cox,1}, T, 'RobustOpts', 'fair', ...
        'ResponseVar', 5, 'Exclude', outl);
    
    % Erstellen des neuen nichtlinearen-Modells ohne die outlier
    Models{3,1} = NonLinearModel.fit(ds{Box_Cox,1}, modelfun, beta0, ...
        'Exclude', outl);
    
    % Print R^2 value of used model
    fprintf('\n\nR^2adj: %1.3f\n\n', Models{ModelChoice,1}.Rsquared.Adjusted)
    
    % Print anzahl outlier
    fprintf('Outlier: %d\n\n', length(outl))
    
    figure(1);
    plotDiagnostics(Models{ModelChoice,1},'cookd') %Cooks Distance
    grid on
    figure(2);
    plotResiduals(Models{ModelChoice,1},'probability') %probability plot
    grid on
    
    bad = input('Nochmal outlier entfernen? Ja(j) oder Nein(n) ', 's');
end
close all

% Outlier entfernen mittels RANSAC nur für Nichtlineares Model
if (~exist('outl', 'var') && (ModelChoice == 3 || ModelChoice == 2))
    rans = input('Sollen outlier mit RANSAC entfernt werden? Ja(j) oder Nein(n) ', 's');

    if (strcmp(rans, 'j'))
        [outl, ~] = RANSAC(14, 10000, 3, 8, ds{Box_Cox}, beta0, T);
        
        % Print anzahl outlier
        fprintf('Outlier: %d\n\n', length(outl))
    end
end

% Anpassen der datasets und Neuberechnung der Modelle
if (exist('outl', 'var'))
    % Löschen der erkannten outlier aus dem Datensatz
    ds{1,1}(outl,:) = [];
    ds{2,1}(outl,:) = [];
    
    % Erstellen eines linearen Regressionsmodells
    Models{1,1} = LinearModel.stepwise(ds{Box_Cox,1}, 'quadratic', ...
        'Lower', 'constant', 'Upper', 'quadratic', 'Verbose', 0);
    
    % Erstellen eines linearen Regressionsmodells
    Models{2,1} = LinearModel.fit(ds{Box_Cox,1}, T, 'RobustOpts', 'fair', ...
        'ResponseVar', 5);
    
    % create non-linear models
    Models{3,1} = NonLinearModel.fit(ds{Box_Cox,1}, modelfun, beta0);
    
    % Print R^2 value of used model
    fprintf('\n\nR^2adj: %1.3f\n\n', Models{ModelChoice,1}.Rsquared.Adjusted)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Diagnose Plots + Vorbereitung auf Plots
% Definieren der benötigten Variablen für confidence intervals Plot
PlotVars.X = linspace(1,length(Models{ModelChoice,1}.Coefficients.Estimate),...
    length(Models{ModelChoice,1}.Coefficients.Estimate))';
PlotVars.Y = Models{ModelChoice,1}.Coefficients.Estimate;
PlotVars.interval = coefCI(Models{ModelChoice,1});
% Plot confidence intervals for the coefficients of the model
figure;
bar(Models{ModelChoice,1}.Coefficients.Estimate)
hold on
errorbar(PlotVars.X, PlotVars.Y, PlotVars.interval(:,1), ...
    PlotVars.interval(:,2),'xr')
xlabel('Koeffizienten')
set(gca,'XTicklabel', Models{ModelChoice,1}.CoefficientNames)
hold off


% Display used formula
disp(Models{ModelChoice,1}.Formula)

% Mittelwerte aus den Bereichen aller Kovariablen
PlotVars.mittel_a = mean([0.2, 0.75]);
PlotVars.mittel_a = repmat(PlotVars.mittel_a, [length(ds{Box_Cox,1}) 1]); %auf vollständige Größe von VP erweitern

PlotVars.mittel_f = mean([200, 2000]);
PlotVars.mittel_f = repmat(PlotVars.mittel_f, [length(ds{Box_Cox,1}) 1]); %auf vollständige Größe von VP erweitern

PlotVars.mittel_v = mean([16.67, 150]);
PlotVars.mittel_v = repmat(PlotVars.mittel_v, [length(ds{Box_Cox,1}) 1]); %auf vollständige Größe von VP erweitern

PlotVars.mittel_p = mean([1000, 3000]);
PlotVars.mittel_p = repmat(PlotVars.mittel_p, [length(ds{Box_Cox,1}) 1]); %auf vollständige Größe von VP erweitern

% Ändern des Namens der Zielvariable für schöne Anzeige in Plot
remain = ziel;
label = char;
while true
    [str, remain] = strtok(remain, '_'); %#ok<STTOK>
    if isempty(str)
        break;
    end
    label = strcat(label, '-', str, '-');
end

%% 2D Scatter Plots
% Warnung falls Box-Cox verwendet wurde
if (Box_Cox == 1)
    warning('Die Konfidenzintervalle sind nicht verlässlich!!!')
end
%% Amplitude
PlotVars.pcolor = linspace(1000,3000,length(ds{Box_Cox,1})); %Leistung wird farbkodiert
% vfsize = linspace(min(VP(:,3)./VP(:,2)),...
%     max(VP(:,3)./VP(:,2)),length(ds{Box_Cox,1})); %Verhältnis v/F wird durch die Größe der "dots" kodiert
% vfsize = vfsize*50; %Größe der "dots" wird angepasst

% Erstelle Variablen für den predictions plot
PlotVars.xscale = linspace(0.2,0.75,length(ds{Box_Cox,1})); %x-Achse skalieren
PlotVars.Xnew = [linspace(0.2,0.75,length(ds{Box_Cox,1}))', ...
    PlotVars.mittel_f, PlotVars.mittel_v, PlotVars.mittel_p]; %Auswertung der prediction an Xnew
[PlotVars.ypred, PlotVars.yci] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew, 'alpha', confi); %prediction

% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred = (lambda*(1/lambda + PlotVars.ypred)).^(1/lambda);
    PlotVars.yci = (lambda*(1/lambda + PlotVars.yci)).^(1/lambda);
end

figure;
scatter(ds{2,1}.a, double(ds{2,1}(:,5)), [], PlotVars.pcolor, 'filled') %scatter plot der gemessenen Daten
xlabel('Amplitude')
xlim([0.15, 0.8])
ylabel(label)
cb = colorbar;
zlab = get(cb,'ylabel');
set(zlab,'String','Leistung');
hold on
errorbar(PlotVars.xscale, PlotVars.ypred, PlotVars.yci(:,1), PlotVars.yci(:,2), 'xk') %plot der prediction mit confidenz interval
hold off

%% Frequenz
PlotVars.pcolor = linspace(1000,3000,length(ds{Box_Cox,1})); %Leistung wird farbkodiert
% avsize = linspace(min(VP(:,1)./VP(:,3)),...
%     max(VP(:,1)./VP(:,3)),length(ds{Box_Cox,1})); %Verhältnis v/F wird durch die Größe der "dots" kodiert
% avsize = avsize*900; %Größe der "dots" wird angepasst

PlotVars.xscale = linspace(200,2000,length(ds{Box_Cox,1}));
PlotVars.Xnew = [PlotVars.mittel_a, linspace(200,2000,length(ds{Box_Cox,1}))', ...
    PlotVars.mittel_v, PlotVars.mittel_p];
[PlotVars.ypred, PlotVars.yci] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew, 'alpha', confi);

% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred = (lambda*(1/lambda + PlotVars.ypred)).^(1/lambda);
    PlotVars.yci = (lambda*(1/lambda + PlotVars.yci)).^(1/lambda);
end

figure;
scatter(ds{2,1}.f, double(ds{2,1}(:,5)), [], PlotVars.pcolor, 'filled')
xlabel('Frequenz')
ylabel(label)
xlim([120 2050])
cb = colorbar;
zlab = get(cb,'ylabel');
set(zlab,'String','Leistung');
hold on
errorbar(PlotVars.xscale, PlotVars.ypred, PlotVars.yci(:,1), PlotVars.yci(:,2), 'xk') %plot der prediction mit confidenz interval
hold off

%% Vorschub
PlotVars.pcolor = linspace(1000,3000,length(ds{Box_Cox,1})); %Leistung wird farbkodiert
% afsize = linspace(min(VP(:,1)./VP(:,2)),...
%     max(VP(:,1)./VP(:,2)),length(ds{Box_Cox,1})); %Verhältnis v/F wird durch die Größe der "dots" kodiert
% afsize = afsize*9000; %Größe der "dots" wird angepasst

PlotVars.xscale = linspace(16.67,150,length(ds{Box_Cox,1}));
PlotVars.Xnew = [PlotVars.mittel_a, PlotVars.mittel_f, ...
    linspace(16.67,150,length(ds{Box_Cox,1}))', PlotVars.mittel_p];
[PlotVars.ypred, PlotVars.yci] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew, 'alpha', confi);

% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred = (lambda*(1/lambda + PlotVars.ypred)).^(1/lambda);
    PlotVars.yci = (lambda*(1/lambda + PlotVars.yci)).^(1/lambda);
end

figure;
scatter(ds{2,1}.v, double(ds{2,1}(:,5)), [], PlotVars.pcolor, 'filled')
xlabel('Vorschub')
ylabel(label)
xlim([11 155])
cb = colorbar;
zlab = get(cb,'ylabel');
set(zlab,'String','Leistung');
hold on
errorbar(PlotVars.xscale, PlotVars.ypred, PlotVars.yci(:,1), PlotVars.yci(:,2), 'xk') %plot der prediction mit confidenz interval
hold off

%% Leistung
PlotVars.acolor = linspace(0.2, 0.75, length(ds{Box_Cox,1})); %Amplitude wird farbkodiert
%vfsize = linspace(min(VP(:,3)./VP(:,2)),...
%max(VP(:,3)./VP(:,2)),length(ds{Box_Cox,1})); %Verhältnis v/F wird durch die Größe der "dots" kodiert
%vfsize = vfsize*50; %Größe der "dots" wird angepasst

PlotVars.xscale = linspace(1000,3000,length(ds{Box_Cox,1}));
PlotVars.Xnew = [PlotVars.mittel_a, PlotVars.mittel_f, PlotVars.mittel_v, ...
    linspace(1000,3000,length(ds{Box_Cox,1}))'];
[PlotVars.ypred, PlotVars.yci] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew, 'alpha', confi);

% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred = (lambda*(1/lambda + PlotVars.ypred)).^(1/lambda);
    PlotVars.yci = (lambda*(1/lambda + PlotVars.yci)).^(1/lambda);
end

figure;
scatter(ds{2,1}.p, double(ds{2,1}(:,5)), [], PlotVars.acolor, 'filled')
xlabel('Leistung')
ylabel(label)
xlim([900 3100])
cb = colorbar;
zlab = get(cb,'ylabel');
set(zlab,'String','Amplitude');
hold on
errorbar(PlotVars.xscale, PlotVars.ypred, PlotVars.yci(:,1), PlotVars.yci(:,2), 'xk') %plot der prediction mit confidenz interval
hold off

%% 3D Scatter Plots
%% Amplitude/Frequenz
PlotVars.pcolor = linspace(1000,3000,length(ds{Box_Cox,1})); %Leistung wird farbkodiert
%vsize = linspace(16.67, 150, length(ds{Box_Cox,1})); %Vorschub wird durch Größe kodiert
%vsize = vsize/2; %Größe der "dots" wird angepasst

PlotVars.xscale = linspace(0.2,0.75,length(ds{Box_Cox,1}))';
PlotVars.yscale = linspace(200,2000,length(ds{Box_Cox,1}))';

PlotVars.predictsZ = zeros(length(PlotVars.yscale), length(PlotVars.yscale));
for i = 1:length(PlotVars.xscale)
	PlotVars.predictsZ(:,i) = predict(Models{ModelChoice,1}, ...
        [repmat(PlotVars.xscale(i), [length(PlotVars.yscale) 1]), ...
        PlotVars.yscale, PlotVars.mittel_v, PlotVars.mittel_p]);
end

% Rücktrafo Box-Cox
if (Box_Cox == 1)
   PlotVars.predictsZ = arrayfun(box_back, PlotVars.predictsZ);
end

figure;
scatter3(ds{2,1}.a, ds{2,1}.f, double(ds{2,1}(:,5)), [], PlotVars.pcolor, 'filled')
set(gca, 'YDir', 'reverse')
xlabel('Amplitude')
ylabel('Frequenz')
zlabel(label)
cb = colorbar;
zlab = get(cb,'ylabel');
set(zlab,'String','Leistung');
hold on
mesh(PlotVars.xscale, PlotVars.yscale, PlotVars.predictsZ)
set(gca, 'YDir', 'reverse')
xlim([0.1, 0.8])
ylim([150, 2100])
hold off

%% Amplitude/Vorschub
PlotVars.pcolor = linspace(1000,3000,length(ds{Box_Cox,1})); %Leistung wird farbkodiert
%fsize = linspace(200, 2000, length(ds{Box_Cox,1})); %Frequenz wird durch Größe kodiert
%fsize = fsize/30; %Größe der "dots" wird angepasst

PlotVars.xscale = linspace(0.2,0.75,length(ds{Box_Cox,1}))';
PlotVars.yscale = linspace(16.67,150,length(ds{Box_Cox,1}))';

PlotVars.predictsZ = zeros(length(PlotVars.yscale), length(PlotVars.yscale));
for i = 1:length(PlotVars.xscale)
	PlotVars.predictsZ(:,i) = predict(Models{ModelChoice,1}, ... 
       [repmat(PlotVars.xscale(i), [length(PlotVars.yscale) 1]), ... 
       PlotVars.mittel_f, PlotVars.yscale, PlotVars.mittel_p]);
end

% Rücktrafo Box-Cox
if (Box_Cox == 1)
   PlotVars.predictsZ = arrayfun(box_back, PlotVars.predictsZ);
end

figure;
scatter3(ds{2,1}.a, ds{2,1}.v, double(ds{2,1}(:,5)), [], PlotVars.pcolor, 'filled')
set(gca, 'YDir', 'reverse')
xlabel('Amplitude')
ylabel('Vorschub')
zlabel(label)
cb = colorbar;
zlab = get(cb,'ylabel');
set(zlab,'String','Leistung');
hold on
mesh(PlotVars.xscale, PlotVars.yscale, PlotVars.predictsZ)
set(gca, 'YDir', 'reverse')
xlim([0.15, 0.8])
ylim([11, 155])
hold off

%% Amplitude/Leistung
PlotVars.vcolor = linspace(16.67, 150, length(ds{Box_Cox,1})); %Vorschub wird farbkodiert
%fsize = linspace(200, 2000, length(ds{Box_Cox,1})); %Frequenz wird durch Größe kodiert
%fsize = fsize/30; %Größe der "dots" wird angepasst

PlotVars.xscale = linspace(0.2,0.75,length(ds{Box_Cox,1}))';
PlotVars.yscale = linspace(1000,3000,length(ds{Box_Cox,1}))';

PlotVars.predictsZ = zeros(length(PlotVars.yscale), length(PlotVars.yscale));
for i = 1:length(PlotVars.xscale)
	PlotVars.predictsZ(:,i) = predict(Models{ModelChoice,1}, ...
        [repmat(PlotVars.xscale(i), [length(PlotVars.yscale) 1]), ...
        PlotVars.mittel_f, PlotVars.mittel_v, PlotVars.yscale]);
end

% Rücktrafo Box-Cox
if (Box_Cox == 1)
   PlotVars.predictsZ = arrayfun(box_back, PlotVars.predictsZ);
end

figure;
scatter3(ds{2,1}.a, ds{2,1}.p, double(ds{2,1}(:,5)), [], PlotVars.vcolor, 'filled')
set(gca, 'YDir', 'reverse')
xlabel('Amplitude')
ylabel('Leistung')
zlabel(label)
cb = colorbar;
zlab = get(cb,'ylabel');
set(zlab,'String','Vorschub');
hold on
mesh(PlotVars.xscale, PlotVars.yscale, PlotVars.predictsZ)
set(gca, 'YDir', 'reverse')
xlim([0.15, 0.8])
ylim([900, 3100])
hold off

%% Frequenz/Vorschub
PlotVars.pcolor = linspace(1000,3000,length(ds{Box_Cox,1})); %Leistung wird farbkodiert
%asize = linspace(0.2, 0.75, length(ds{Box_Cox,1})); %Amplitude wird durch Größe kodiert
%asize = asize*80; %Größe der "dots" wird angepasst

PlotVars.xscale = linspace(200,2000,length(ds{Box_Cox,1}))';
PlotVars.yscale = linspace(16.67,150,length(ds{Box_Cox,1}))';

PlotVars.predictsZ = zeros(length(PlotVars.yscale), length(PlotVars.yscale));
for i = 1:length(PlotVars.xscale)
	PlotVars.predictsZ(:,i) = predict(Models{ModelChoice,1}, ...
        [PlotVars.mittel_a, repmat(PlotVars.xscale(i), ...
        [length(PlotVars.yscale) 1]), PlotVars.yscale, PlotVars.mittel_p]);
end

% Rücktrafo Box-Cox
if (Box_Cox == 1)
   PlotVars.predictsZ = arrayfun(box_back, PlotVars.predictsZ);
end

figure;
scatter3(ds{2,1}.f, ds{2,1}.v, double(ds{2,1}(:,5)), [], PlotVars.pcolor, 'filled')
set(gca, 'YDir', 'reverse')
xlabel('Frequenz')
ylabel('Vorschub')
zlabel(label)
cb = colorbar;
zlab = get(cb,'ylabel');
set(zlab,'String','Leistung');
hold on
mesh(PlotVars.xscale, PlotVars.yscale, PlotVars.predictsZ)
set(gca, 'YDir', 'reverse')
xlim([120, 2100])
ylim([12, 160])
hold off

%% Frequenz/Leistung
PlotVars.vcolor = linspace(16.67, 150, length(ds{Box_Cox,1})); %Vorschub wird farbkodiert
%asize = linspace(0.2, 0.75, length(ds{Box_Cox,1})); %Amplitude wird durch Größe kodiert
%asize = asize*80; %Größe der "dots" wird angepasst

PlotVars.xscale = linspace(200,2000,length(ds{Box_Cox,1}))';
PlotVars.yscale = linspace(1000,3000,length(ds{Box_Cox,1}))';

PlotVars.predictsZ = zeros(length(PlotVars.yscale), length(PlotVars.yscale));
for i = 1:length(PlotVars.xscale)
	PlotVars.predictsZ(:,i) = predict(Models{ModelChoice,1}, ...
        [PlotVars.mittel_a, repmat(PlotVars.xscale(i), ...
        [length(PlotVars.yscale) 1]), PlotVars.mittel_v, PlotVars.yscale]);
end

% Rücktrafo Box-Cox
if (Box_Cox == 1)
   PlotVars.predictsZ = arrayfun(box_back, PlotVars.predictsZ);
end

figure;
scatter3(ds{2,1}.f, ds{2,1}.p, double(ds{2,1}(:,5)), [], PlotVars.vcolor, 'filled')
set(gca, 'YDir', 'reverse')
xlabel('Frequenz')
ylabel('Leistung')
zlabel(label)
cb = colorbar;
zlab = get(cb,'ylabel');
set(zlab,'String','Vorschub');
hold on
mesh(PlotVars.xscale, PlotVars.yscale, PlotVars.predictsZ)
set(gca, 'YDir', 'reverse')
xlim([120, 2100])
ylim([900, 3100])
hold off

%% Vorschub/Leistung
PlotVars.fcolor = linspace(200, 2000, length(ds{Box_Cox,1})); %Frequenz wird farbkodiert
%asize = linspace(0.2, 0.75, length(ds{Box_Cox,1})); %Amplitude wird durch Größe kodiert
%asize = asize*80; %Größe der "dots" wird angepasst

PlotVars.xscale = linspace(16.67,150,length(ds{Box_Cox,1}))';
PlotVars.yscale = linspace(1000,3000,length(ds{Box_Cox,1}))';

PlotVars.predictsZ = zeros(length(PlotVars.yscale), length(PlotVars.yscale));
for i = 1:length(PlotVars.xscale)
	PlotVars.predictsZ(:,i) = predict(Models{ModelChoice,1}, ...
        [PlotVars.mittel_a, PlotVars.mittel_f, repmat(PlotVars.xscale(i), ...
        [length(PlotVars.yscale) 1]), PlotVars.yscale]);
end

% Rücktrafo Box-Cox
if (Box_Cox == 1)
   PlotVars.predictsZ = arrayfun(box_back, PlotVars.predictsZ);
end

figure;
scatter3(ds{2,1}.v, ds{2,1}.p, double(ds{2,1}(:,5)), [], PlotVars.fcolor, 'filled')
set(gca, 'YDir', 'reverse')
xlabel('Vorschub')
ylabel('Leistung')
zlabel(label)
cb = colorbar;
zlab = get(cb,'ylabel');
set(zlab,'String','Frequenz');
hold on
mesh(PlotVars.xscale, PlotVars.yscale, PlotVars.predictsZ)
set(gca, 'YDir', 'reverse')
xlim([12, 160])
ylim([900, 3100])
hold off

%% Zusammenfassung des Modells
Information = cell(17,1);
if (ModelChoice == 3)
    Information{1,1} = strcat('Modellgleichung: ', Models{ModelChoice,1}.Formula.Expression);
else
    Information{1,1} = strcat('Modellgleichung: ', Models{ModelChoice,1}.Formula.LinearPredictor);
end
Information{2,1} = '';
Information{3,1} = 'Werkstoff: Edelstahl';
Information{4,1} = '';
Information{5,1} = strcat('R^2Adj: ', num2str(Models{ModelChoice,1}.Rsquared.Adjusted, 3));
Information{6,1} = '';
Information{7,1} = strcat('R^2: ', num2str(Models{ModelChoice,1}.Rsquared.Ordinary, 3));
Information{8,1} = '';
Information{9,1} = strcat('AIC: ', num2str(Models{ModelChoice,1}.ModelCriterion.AIC, 3));
Information{10,1} = '';
if (exist('outl', 'var'))
    Information{11,1} = strcat('Anzahl Außreiser: ', int2str(length(outl)));
end
Information{12,1} = '';
Information{13,1} = strcat('VP: ', int2str(length(VP)));
Information{14,1} = '';
Information{15,1} = strcat('Filter: ', filter);
Information{14,1} = '';
if (exist('outl', 'var'))
    Information{16,1} = strcat('outl: ', int2str(outl));
end

% Definieren der benötigten Variablen für confidence intervals Plot
PlotVars.X = linspace(1,length(Models{ModelChoice,1}.Coefficients.Estimate),...
    length(Models{ModelChoice,1}.Coefficients.Estimate))';
PlotVars.Y = Models{ModelChoice,1}.Coefficients.Estimate;
PlotVars.interval = coefCI(Models{ModelChoice,1});

% Plot
figure('Name','Zusammenfassung','NumberTitle','off','units',...
    'normalized','outerposition',[0 0 1 1]);
subplot(2,2,3), plotResiduals(Models{ModelChoice,1},'probability')
subplot(2,2,4),
bar(Models{ModelChoice,1}.Coefficients.Estimate)
hold on
errorbar(PlotVars.X,PlotVars.Y,PlotVars.interval(:,1),PlotVars.interval(:,2),'xr')
xlabel('Koeffizienten')
set(gca,'XTicklabel', Models{ModelChoice,1}.CoefficientNames)
hold off
subplot(2,2,[1,2]),
axis off
annotation('textbox', [.38, .85, .1, .1],'String', Information)
set(gca,'ButtonDownFcn','selectmoveresize');
set(gcf,'WindowButtonDownFcn','selectmoveresize');

if (Box_Cox == 2)
    plotSlice(Models{ModelChoice,1})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plots für VÖ

% für f_1 --> v<=20
% für f_2 --> v gesamt

% für v_1 --> f gesamt
% für v_2 --> f>=1500

% für a_1 --> f komplett
% für a_2 --> f<=1449

% für f_1 --> a gesamt
% für f_2 = 1743 --> a<=0.38

if (strcmp(ziel, 'VeloDis_E') || strcmp(ziel, 'VeloDis_B'))
    PlotVars.Umrechnung = 0.1035;
else
    PlotVars.Umrechnung = 1;
end

PlotVars.yLABEL = num2cell(0:PlotVars.Step:PlotVars.IntyScale);
PlotVars.yLABEL = cellfun(@num2str, PlotVars.yLABEL, 'UniformOutput', 0);
PlotVars.yLABEL{1,length(PlotVars.yLABEL)-1} = 'm/s';

PlotVars.xLABEL = num2cell(1000:1000:3000);
PlotVars.xLABEL = cellfun(@num2str, PlotVars.xLABEL, 'UniformOutput', 0);
PlotVars.xLABEL{1,length(PlotVars.xLABEL)-1} = 'W';

% Plot
figure('Name','Interactionplot','NumberTitle','off','units',...
    'centimeters');
subplot(4,4,1), axis off
%%%
% Berechnung der für Plot benötigten Variablen
%%% a:f %%%
PlotVars.a_1 = 0.2;
PlotVars.a_2 = 0.75;
PlotVars.a_1 = repmat(PlotVars.a_1, [length(ds{Box_Cox,1}) 1]);
PlotVars.a_2 = repmat(PlotVars.a_2, [length(ds{Box_Cox,1}) 1]);
PlotVars.xscale = linspace(200,2000,length(ds{Box_Cox,1}));
PlotVars.Xnew_1 = [PlotVars.a_1, linspace(200,2000,length(ds{Box_Cox,1}))', ...
    PlotVars.mittel_v, PlotVars.mittel_p];
PlotVars.Xnew_2 = [PlotVars.a_2, linspace(200,2000,length(ds{Box_Cox,1}))', ...
    PlotVars.mittel_v, PlotVars.mittel_p];
[PlotVars.ypred_1, PlotVars.yci_1] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_1, 'alpha', confi);
[PlotVars.ypred_2, PlotVars.yci_2] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_2, 'alpha', confi);
% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred_1 = (lambda*(1/lambda + PlotVars.ypred_1)).^(1/lambda);
    PlotVars.ypred_2 = (lambda*(1/lambda + PlotVars.ypred_2)).^(1/lambda);
    PlotVars.yci_1 = (lambda*(1/lambda + PlotVars.yci_1)).^(1/lambda);
    PlotVars.yci_2 = (lambda*(1/lambda + PlotVars.yci_2)).^(1/lambda);
end
%%%
% Einschränkung des VP --> Kürzung der Arrays
%%%
cut = find(PlotVars.Xnew_2(:,2) < 1449);
PlotVars.ypred_2 = PlotVars.ypred_2(1:cut(end));
PlotVars.yci_2 = PlotVars.yci_2(1:cut(end),:);
PlotVars.xscale_a = PlotVars.xscale(1:cut(end));
%%%
subplot(4,4,2), plot(PlotVars.xscale, PlotVars.ypred_1*PlotVars.Umrechnung, 'Color', [1 0.5 0])
hold on
plot(PlotVars.xscale_a, PlotVars.ypred_2*PlotVars.Umrechnung, 'b')
plot(PlotVars.xscale, PlotVars.yci_1*PlotVars.Umrechnung, '--', 'Color', [1 0.5 0])
plot(PlotVars.xscale_a, PlotVars.yci_2*PlotVars.Umrechnung, '--b')
set(gca,'units','centimeters','Position',[4.5, 11.5, 3, 3], 'Fontsize', FontSize)
h = legend('a = 0.2 mm', 'a = 0.75 mm');
set(h, 'units', 'centimeters', 'position', ...
    [0.49 ,12.49767364583333 ,2.8, 1.0046527083333332], ...
    'Fontsize', FontSize)
hold off
g = ylabel('Zielgröße', 'Fontsize', FontSize+2);
set(g, 'Units', 'centimeters', 'Position', [-0.4, 1.4805408333333332, 0])
set(gca, 'XTickLabel', {})
set(gca, 'YTick', 0:PlotVars.Step:PlotVars.IntyScale)
set(gca, 'YTickLabel', PlotVars.yLABEL)
ylim([0 PlotVars.IntyScale])
xlim([200 2000])
%%%
% Berechnung der für Plot benötigten Variablen
%%% a:v %%%
PlotVars.xscale = linspace(16.67,150,length(ds{Box_Cox,1}));
PlotVars.Xnew_1 = [PlotVars.a_1, PlotVars.mittel_f, ...
    linspace(16.67,150,length(ds{Box_Cox,1}))', PlotVars.mittel_p];
PlotVars.Xnew_2 = [PlotVars.a_2, PlotVars.mittel_f, ...
    linspace(16.67,150,length(ds{Box_Cox,1}))', PlotVars.mittel_p];
[PlotVars.ypred_1, PlotVars.yci_1] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_1, 'alpha', confi);
[PlotVars.ypred_2, PlotVars.yci_2] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_2, 'alpha', confi);
% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred_1 = (lambda*(1/lambda + PlotVars.ypred_1)).^(1/lambda);
    PlotVars.ypred_2 = (lambda*(1/lambda + PlotVars.ypred_2)).^(1/lambda);
    PlotVars.yci_1 = (lambda*(1/lambda + PlotVars.yci_1)).^(1/lambda);
    PlotVars.yci_2 = (lambda*(1/lambda + PlotVars.yci_2)).^(1/lambda);
end
%%%
subplot(4,4,3), plot(PlotVars.xscale, PlotVars.ypred_1*PlotVars.Umrechnung, 'Color', [1 0.5 0])
hold on
plot(PlotVars.xscale, PlotVars.ypred_2*PlotVars.Umrechnung, 'b')
plot(PlotVars.xscale, PlotVars.yci_1*PlotVars.Umrechnung, '--', 'Color', [1 0.5 0])
plot(PlotVars.xscale, PlotVars.yci_2*PlotVars.Umrechnung, '--b')
set(gca,'units','centimeters','Position',[8, 11.5, 3, 3], 'Fontsize', FontSize)
set(gca, 'XTickLabel', {})
set(gca, 'YTick', 0:PlotVars.Step:PlotVars.IntyScale)
set(gca, 'YTickLabel', {})
ylim([0 PlotVars.IntyScale])
xlim([16.67 150])
%%%
% Berechnung der für Plot benötigten Variablen
%%% a:p %%%
PlotVars.xscale = linspace(1000,3000,length(ds{Box_Cox,1}));
PlotVars.Xnew_1 = [PlotVars.a_1, PlotVars.mittel_f, PlotVars.mittel_v, ...
    linspace(1000,3000,length(ds{Box_Cox,1}))'];
PlotVars.Xnew_2 = [PlotVars.a_2, PlotVars.mittel_f, PlotVars.mittel_v, ...
    linspace(1000,3000,length(ds{Box_Cox,1}))'];
[PlotVars.ypred_1, PlotVars.yci_1] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_1, 'alpha', confi);
[PlotVars.ypred_2, PlotVars.yci_2] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_2, 'alpha', confi);
% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred_1 = (lambda*(1/lambda + PlotVars.ypred_1)).^(1/lambda);
    PlotVars.ypred_2 = (lambda*(1/lambda + PlotVars.ypred_2)).^(1/lambda);
    PlotVars.yci_1 = (lambda*(1/lambda + PlotVars.yci_1)).^(1/lambda);
    PlotVars.yci_2 = (lambda*(1/lambda + PlotVars.yci_2)).^(1/lambda);
end
%%%
subplot(4,4,4), plot(PlotVars.xscale, PlotVars.ypred_1*PlotVars.Umrechnung, 'Color', [1 0.5 0])
hold on
plot(PlotVars.xscale, PlotVars.ypred_2*PlotVars.Umrechnung, 'b')
plot(PlotVars.xscale, PlotVars.yci_1*PlotVars.Umrechnung, '--', 'Color', [1 0.5 0])
plot(PlotVars.xscale, PlotVars.yci_2*PlotVars.Umrechnung, '--b')
set(gca,'units','centimeters','Position',[11.5, 11.5, 3, 3], 'Fontsize', FontSize)
set(gca, 'XTickLabel', {})
set(gca, 'YTick', 0:PlotVars.Step:PlotVars.IntyScale)
set(gca, 'YTickLabel', {})
ylim([0 PlotVars.IntyScale])
xlim([1000 3000])
%%%
% Berechnung der für Plot benötigten Variablen
%%% f:a %%%
PlotVars.f_1 = 200;
PlotVars.f_2 = 1743;
PlotVars.f_1 = repmat(PlotVars.f_1, [length(ds{Box_Cox,1}) 1]);
PlotVars.f_2 = repmat(PlotVars.f_2, [length(ds{Box_Cox,1}) 1]);
PlotVars.xscale = linspace(0.2,0.75,length(ds{Box_Cox,1}));
PlotVars.Xnew_1 = [linspace(0.2,0.75,length(ds{Box_Cox,1}))', ...
    PlotVars.f_1, PlotVars.mittel_v, PlotVars.mittel_p];
PlotVars.Xnew_2 = [linspace(0.2,0.75,length(ds{Box_Cox,1}))', ...
    PlotVars.f_2, PlotVars.mittel_v, PlotVars.mittel_p];
[PlotVars.ypred_1, PlotVars.yci_1] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_1, 'alpha', confi);
[PlotVars.ypred_2, PlotVars.yci_2] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_2, 'alpha', confi);
% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred_1 = (lambda*(1/lambda + PlotVars.ypred_1)).^(1/lambda);
    PlotVars.ypred_2 = (lambda*(1/lambda + PlotVars.ypred_2)).^(1/lambda);
    PlotVars.yci_1 = (lambda*(1/lambda + PlotVars.yci_1)).^(1/lambda);
    PlotVars.yci_2 = (lambda*(1/lambda + PlotVars.yci_2)).^(1/lambda);
end
%%%
% Einschränkung des VP --> Kürzung der Arrays
%%%
cut = find(PlotVars.Xnew_2(:,1) < 0.38);
PlotVars.ypred_2 = PlotVars.ypred_2(1:cut(end));
PlotVars.yci_2 = PlotVars.yci_2(1:cut(end),:);
PlotVars.xscale_a = PlotVars.xscale(1:cut(end));
%%%
subplot(4,4,5), plot(PlotVars.xscale, PlotVars.ypred_1*PlotVars.Umrechnung, 'Color', [1 0.5 0])
hold on
plot(PlotVars.xscale_a, PlotVars.ypred_2*PlotVars.Umrechnung, 'b')
plot(PlotVars.xscale, PlotVars.yci_1*PlotVars.Umrechnung, '--', 'Color', [1 0.5 0])
plot(PlotVars.xscale_a, PlotVars.yci_2*PlotVars.Umrechnung, '--b')
set(gca,'units','centimeters','Position',[1, 8, 3, 3],'Fontsize', FontSize)
h = legend('f = 200 Hz', 'f = 1743 Hz');
set(h, 'Units', 'centimeters', 'Position', ...
    [4.3, 8.99767364583333, 3.4, 1.0046527083333332], ...
    'Fontsize', FontSize)
g = ylabel('Zielgröße', 'Fontsize', FontSize+2);
set(g, 'Units', 'centimeters', 'Position', [-0.4, 1.4805408333333332, 0])
set(gca, 'XTickLabel', {})
set(gca, 'YTickLabel', PlotVars.yLABEL)
set(gca, 'YTick', 0:PlotVars.Step:PlotVars.IntyScale)
ylim([0 PlotVars.IntyScale])
xlim([0.2 0.75])
%%%
subplot(4,4,6), axis off
%%%
% Berechnung der für Plot benötigten Variablen
%%% f:v %%%
PlotVars.xscale = linspace(16.67,150,length(ds{Box_Cox,1}));
PlotVars.Xnew_1 = [PlotVars.mittel_a, PlotVars.f_1, ...
    linspace(16.67,150,length(ds{Box_Cox,1}))', PlotVars.mittel_p];
PlotVars.Xnew_2 = [PlotVars.mittel_a, PlotVars.f_2, ...
    linspace(16.67,150,length(ds{Box_Cox,1}))', PlotVars.mittel_p];
[PlotVars.ypred_1, PlotVars.yci_1] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_1, 'alpha', confi);
[PlotVars.ypred_2, PlotVars.yci_2] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_2, 'alpha', confi);
% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred_1 = (lambda*(1/lambda + PlotVars.ypred_1)).^(1/lambda);
    PlotVars.ypred_2 = (lambda*(1/lambda + PlotVars.ypred_2)).^(1/lambda);
    PlotVars.yci_1 = (lambda*(1/lambda + PlotVars.yci_1)).^(1/lambda);
    PlotVars.yci_2 = (lambda*(1/lambda + PlotVars.yci_2)).^(1/lambda);
end
%%%
% Einschränkung des VP --> Kürzung der Arrays
%%%
cut = find(PlotVars.Xnew_1(:,3) < 20);
PlotVars.ypred_1 = PlotVars.ypred_1(1:cut(end));
PlotVars.yci_1 = PlotVars.yci_1(1:cut(end),:);
PlotVars.xscale_v = PlotVars.xscale(1:cut(end));
%%%
subplot(4,4,7), plot(PlotVars.xscale_v, PlotVars.ypred_1*PlotVars.Umrechnung, 'Color', [1 0.5 0])
hold on
plot(PlotVars.xscale, PlotVars.ypred_2*PlotVars.Umrechnung, 'b')
plot(PlotVars.xscale_v, PlotVars.yci_1*PlotVars.Umrechnung, '--', 'Color', [1 0.5 0])
plot(PlotVars.xscale, PlotVars.yci_2*PlotVars.Umrechnung, '--b')
set(gca,'units','centimeters','Position',[8, 8, 3, 3], 'Fontsize', FontSize)
set(gca, 'XTickLabel', {})
set(gca, 'YTickLabel', {})
set(gca, 'YTick', 0:PlotVars.Step:PlotVars.IntyScale)
ylim([0 PlotVars.IntyScale])
xlim([16.67 150])
%%%
% Berechnung der für Plot benötigten Variablen
%%% f:p %%%
PlotVars.xscale = linspace(1000,3000,length(ds{Box_Cox,1}));
PlotVars.Xnew_1 = [PlotVars.mittel_a, PlotVars.f_1, PlotVars.mittel_v, ...
    linspace(1000,3000,length(ds{Box_Cox,1}))'];
PlotVars.Xnew_2 = [PlotVars.mittel_a, PlotVars.f_2, PlotVars.mittel_v, ...
    linspace(1000,3000,length(ds{Box_Cox,1}))'];
[PlotVars.ypred_1, PlotVars.yci_1] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_1, 'alpha', confi);
[PlotVars.ypred_2, PlotVars.yci_2] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_2, 'alpha', confi);
% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred_1 = (lambda*(1/lambda + PlotVars.ypred_1)).^(1/lambda);
    PlotVars.ypred_2 = (lambda*(1/lambda + PlotVars.ypred_2)).^(1/lambda);
    PlotVars.yci_1 = (lambda*(1/lambda + PlotVars.yci_1)).^(1/lambda);
    PlotVars.yci_2 = (lambda*(1/lambda + PlotVars.yci_2)).^(1/lambda);
end
%%%
subplot(4,4,8), plot(PlotVars.xscale, PlotVars.ypred_1*PlotVars.Umrechnung, 'Color', [1 0.5 0])
hold on
plot(PlotVars.xscale, PlotVars.ypred_2*PlotVars.Umrechnung, 'b')
plot(PlotVars.xscale, PlotVars.yci_1*PlotVars.Umrechnung, '--', 'Color', [1 0.5 0])
plot(PlotVars.xscale, PlotVars.yci_2*PlotVars.Umrechnung, '--b')
set(gca,'units','centimeters','Position',[11.5, 8, 3, 3], 'Fontsize', FontSize)
set(gca, 'XTickLabel', {})
set(gca, 'YTickLabel', {})
set(gca, 'YTick', 0:PlotVars.Step:PlotVars.IntyScale)
ylim([0 PlotVars.IntyScale])
xlim([1000 3000])
%%%
% Berechnung der für Plot benötigten Variablen
%%% v:a %%%
PlotVars.v_1 = 16.67;
PlotVars.v_2 = 150;
PlotVars.v_1 = repmat(PlotVars.v_1, [length(ds{Box_Cox,1}) 1]);
PlotVars.v_2 = repmat(PlotVars.v_2, [length(ds{Box_Cox,1}) 1]);
PlotVars.xscale = linspace(0.2,0.75,length(ds{Box_Cox,1}));
PlotVars.Xnew_1 = [linspace(0.2,0.75,length(ds{Box_Cox,1}))', ...
    PlotVars.mittel_f, PlotVars.v_1, PlotVars.mittel_p];
PlotVars.Xnew_2 = [linspace(0.2,0.75,length(ds{Box_Cox,1}))', ...
    PlotVars.mittel_f, PlotVars.v_2, PlotVars.mittel_p];
[PlotVars.ypred_1, PlotVars.yci_1] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_1, 'alpha', confi);
[PlotVars.ypred_2, PlotVars.yci_2] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_2, 'alpha', confi);
% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred_1 = (lambda*(1/lambda + PlotVars.ypred_1)).^(1/lambda);
    PlotVars.ypred_2 = (lambda*(1/lambda + PlotVars.ypred_2)).^(1/lambda);
    PlotVars.yci_1 = (lambda*(1/lambda + PlotVars.yci_1)).^(1/lambda);
    PlotVars.yci_2 = (lambda*(1/lambda + PlotVars.yci_2)).^(1/lambda);
end
%%%
subplot(4,4,9), plot(PlotVars.xscale, PlotVars.ypred_1*PlotVars.Umrechnung, 'Color', [1 0.5 0])
hold on
plot(PlotVars.xscale, PlotVars.ypred_2*PlotVars.Umrechnung, 'b')
plot(PlotVars.xscale, PlotVars.yci_1*PlotVars.Umrechnung, '--', 'Color', [1 0.5 0])
plot(PlotVars.xscale, PlotVars.yci_2*PlotVars.Umrechnung, '--b')
set(gca,'units','centimeters','Position',[1, 4.5, 3, 3], 'Fontsize', FontSize)
h = legend('v = 16,67 mm/s', 'v = 150 mm/s');
set(h, 'Units', 'centimeters', 'Position', ...
    [7.768295989583338, 5.497673645833332, 3.4, 1.0046527083333332], ...
    'Fontsize', FontSize)
set(gca, 'XTickLabel', {})
set(gca, 'YTickLabel', PlotVars.yLABEL)
set(gca, 'YTick', 0:PlotVars.Step:PlotVars.IntyScale)
g = ylabel('Zielgröße', 'Fontsize', FontSize+2);
set(g, 'Units', 'centimeters', 'Position', [-0.4, 1.4805408333333332, 0])
ylim([0 PlotVars.IntyScale])
xlim([0.2 0.75])
%%%
% Berechnung der für Plot benötigten Variablen
%%% v:f %%%
PlotVars.xscale = linspace(200,2000,length(ds{Box_Cox,1}));
PlotVars.Xnew_1 = [PlotVars.mittel_a, linspace(200,2000,length(ds{Box_Cox,1}))', ...
    PlotVars.v_1, PlotVars.mittel_p];
PlotVars.Xnew_2 = [PlotVars.mittel_a, linspace(200,2000,length(ds{Box_Cox,1}))', ...
    PlotVars.v_2, PlotVars.mittel_p];
[PlotVars.ypred_1, PlotVars.yci_1] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_1, 'alpha', confi);
[PlotVars.ypred_2, PlotVars.yci_2] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_2, 'alpha', confi);
% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred_1 = (lambda*(1/lambda + PlotVars.ypred_1)).^(1/lambda);
    PlotVars.ypred_2 = (lambda*(1/lambda + PlotVars.ypred_2)).^(1/lambda);
    PlotVars.yci_1 = (lambda*(1/lambda + PlotVars.yci_1)).^(1/lambda);
    PlotVars.yci_2 = (lambda*(1/lambda + PlotVars.yci_2)).^(1/lambda);
end
%%%
% Einschränkung des VP --> Kürzung der Arrays
%%%
cut = find(PlotVars.Xnew_2(:,2) > 1500);
PlotVars.ypred_2 = PlotVars.ypred_2(cut(1):cut(end));
PlotVars.yci_2 = PlotVars.yci_2(cut(1):cut(end),:);
PlotVars.xscale_f = PlotVars.xscale(cut(1):cut(end));
%%%
subplot(4,4,10), plot(PlotVars.xscale, PlotVars.ypred_1*PlotVars.Umrechnung, 'Color', [1 0.5 0])
hold on
plot(PlotVars.xscale_f, PlotVars.ypred_2*PlotVars.Umrechnung, 'b')
plot(PlotVars.xscale, PlotVars.yci_1*PlotVars.Umrechnung, '--', 'Color', [1 0.5 0])
plot(PlotVars.xscale_f, PlotVars.yci_2*PlotVars.Umrechnung, '--b')
set(gca,'units','centimeters','Position',[4.5, 4.5, 3, 3], 'Fontsize', FontSize)
set(gca, 'XTickLabel', {})
set(gca, 'YTickLabel', {})
set(gca, 'YTick', 0:PlotVars.Step:PlotVars.IntyScale)
ylim([0 PlotVars.IntyScale])
xlim([200 2000])
%%%
subplot(4,4,11), axis off
%%%
% Berechnung der für Plot benötigten Variablen
%%% v:p %%%
PlotVars.xscale = linspace(1000,3000,length(ds{Box_Cox,1}));
PlotVars.Xnew_1 = [PlotVars.mittel_a, PlotVars.mittel_f, PlotVars.v_1, ...
    linspace(1000,3000,length(ds{Box_Cox,1}))'];
PlotVars.Xnew_2 = [PlotVars.mittel_a, PlotVars.mittel_f, PlotVars.v_2, ...
    linspace(1000,3000,length(ds{Box_Cox,1}))'];
[PlotVars.ypred_1, PlotVars.yci_1] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_1, 'alpha', confi);
[PlotVars.ypred_2, PlotVars.yci_2] = predict(Models{ModelChoice,1}, ...
    PlotVars.Xnew_2, 'alpha', confi);
% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred_1 = (lambda*(1/lambda + PlotVars.ypred_1)).^(1/lambda);
    PlotVars.ypred_2 = (lambda*(1/lambda + PlotVars.ypred_2)).^(1/lambda);
    PlotVars.yci_1 = (lambda*(1/lambda + PlotVars.yci_1)).^(1/lambda);
    PlotVars.yci_2 = (lambda*(1/lambda + PlotVars.yci_2)).^(1/lambda);
end
%%%
subplot(4,4,12), plot(PlotVars.xscale, PlotVars.ypred_1*PlotVars.Umrechnung, 'Color', [1 0.5 0])
hold on
plot(PlotVars.xscale, PlotVars.ypred_2*PlotVars.Umrechnung, 'b')
plot(PlotVars.xscale, PlotVars.yci_1*PlotVars.Umrechnung, '--', 'Color', [1 0.5 0])
plot(PlotVars.xscale, PlotVars.yci_2*PlotVars.Umrechnung, '--b')
set(gca,'units','centimeters','Position',[11.5, 4.5, 3, 3], 'Fontsize', FontSize)
xlabel('p', 'Fontsize', FontSize+2)
set(gca, 'YTickLabel', {})
set(gca, 'XTickLabel', PlotVars.xLABEL)
set(gca, 'YTick', 0:PlotVars.Step:PlotVars.IntyScale)
ylim([0 PlotVars.IntyScale])
xlim([1000 3000])
%%%
% Berechnung der für Plot benötigten Variablen
%%% p:a %%%
PlotVars.p_1 = 1000;
PlotVars.p_2 = 3000;
PlotVars.p_1 = repmat(PlotVars.p_1, [length(ds{Box_Cox,1}) 1]);
PlotVars.p_2 = repmat(PlotVars.p_2, [length(ds{Box_Cox,1}) 1]);
PlotVars.xscale = linspace(0.2,0.75,length(ds{Box_Cox,1}));
PlotVars.Xnew_1 = [linspace(0.2,0.75,length(ds{Box_Cox,1}))', PlotVars.mittel_f, PlotVars.mittel_v, PlotVars.p_1];
PlotVars.Xnew_2 = [linspace(0.2,0.75,length(ds{Box_Cox,1}))', PlotVars.mittel_f, PlotVars.mittel_v, PlotVars.p_2];
[PlotVars.ypred_1, PlotVars.yci_1] = predict(Models{ModelChoice,1}, PlotVars.Xnew_1, 'alpha', confi);
[PlotVars.ypred_2, PlotVars.yci_2] = predict(Models{ModelChoice,1}, PlotVars.Xnew_2, 'alpha', confi);
% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred_1 = (lambda*(1/lambda + PlotVars.ypred_1)).^(1/lambda);
    PlotVars.ypred_2 = (lambda*(1/lambda + PlotVars.ypred_2)).^(1/lambda);
    PlotVars.yci_1 = (lambda*(1/lambda + PlotVars.yci_1)).^(1/lambda);
    PlotVars.yci_2 = (lambda*(1/lambda + PlotVars.yci_2)).^(1/lambda);
end
%%%
PlotVars.xLABEL = num2cell(linspace(0.25,0.75,3));
PlotVars.xLABEL = cellfun(@num2str, PlotVars.xLABEL, 'UniformOutput', 0);
PlotVars.xLABEL{1,length(PlotVars.xLABEL)-1} = 'mm';

subplot(4,4,13), plot(PlotVars.xscale, PlotVars.ypred_1*PlotVars.Umrechnung, 'Color', [1 0.5 0])
hold on
plot(PlotVars.xscale, PlotVars.ypred_2*PlotVars.Umrechnung, 'b')
plot(PlotVars.xscale, PlotVars.yci_1*PlotVars.Umrechnung, '--', 'Color', [1 0.5 0])
plot(PlotVars.xscale, PlotVars.yci_2*PlotVars.Umrechnung, '--b')
set(gca,'units','centimeters','Position',[1, 1, 3, 3], 'Fontsize', FontSize)
h = legend('p = 1000 W', 'p = 3000 W');
set(h, 'Units', 'centimeters', 'Position', ...
    [11.506240052083337, 1.971235416666667, 3.4, 1.0046527083333332], ...
    'Fontsize', FontSize)
xlabel('a', 'Fontsize', FontSize+2)
g = ylabel('Zielgröße', 'Fontsize', FontSize+2);
set(g, 'Units', 'centimeters', 'Position', [-0.4, 1.4805408333333332, 0])
set(gca, 'YTickLabel', PlotVars.yLABEL)
set(gca, 'XTick', linspace(0.25,0.75,3))
set(gca, 'XTickLabel', PlotVars.xLABEL)
set(gca, 'YTick', 0:PlotVars.Step:PlotVars.IntyScale)
ylim([0 PlotVars.IntyScale])
xlim([0.2 0.75])
%%%
% Berechnung der für Plot benötigten Variablen
%%% p:f %%%
PlotVars.xscale = linspace(200,2000,length(ds{Box_Cox,1}));
PlotVars.Xnew_1 = [PlotVars.mittel_a, linspace(200,2000,length(ds{Box_Cox,1}))', PlotVars.mittel_v, PlotVars.p_1];
PlotVars.Xnew_2 = [PlotVars.mittel_a, linspace(200,2000,length(ds{Box_Cox,1}))', PlotVars.mittel_v, PlotVars.p_2];
[PlotVars.ypred_1, PlotVars.yci_1] = predict(Models{ModelChoice,1}, PlotVars.Xnew_1, 'alpha', confi);
[PlotVars.ypred_2, PlotVars.yci_2] = predict(Models{ModelChoice,1}, PlotVars.Xnew_2, 'alpha', confi);
% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred_1 = (lambda*(1/lambda + PlotVars.ypred_1)).^(1/lambda);
    PlotVars.ypred_2 = (lambda*(1/lambda + PlotVars.ypred_2)).^(1/lambda);
    PlotVars.yci_1 = (lambda*(1/lambda + PlotVars.yci_1)).^(1/lambda);
    PlotVars.yci_2 = (lambda*(1/lambda + PlotVars.yci_2)).^(1/lambda);
end
%%%
PlotVars.xLABEL = num2cell(500:500:2000);
PlotVars.xLABEL = cellfun(@num2str, PlotVars.xLABEL, 'UniformOutput', 0);
PlotVars.xLABEL{1,length(PlotVars.xLABEL)-1} = 'Hz';

subplot(4,4,14), plot(PlotVars.xscale, PlotVars.ypred_1*PlotVars.Umrechnung, 'Color', [1 0.5 0])
hold on
plot(PlotVars.xscale, PlotVars.ypred_2*PlotVars.Umrechnung, 'b')
plot(PlotVars.xscale, PlotVars.yci_1*PlotVars.Umrechnung, '--', 'Color', [1 0.5 0])
plot(PlotVars.xscale, PlotVars.yci_2*PlotVars.Umrechnung, '--b')
set(gca,'units','centimeters','Position',[4.5, 1, 3, 3], 'Fontsize', FontSize)
xlabel('f', 'Fontsize', FontSize+2)
set(gca, 'YTickLabel', {})
set(gca, 'XTick', 500:500:2000)
set(gca, 'XTickLabel', PlotVars.xLABEL)
set(gca, 'YTick', 0:PlotVars.Step:PlotVars.IntyScale)
ylim([0 PlotVars.IntyScale])
xlim([200 2000])
% Berechnung der für Plot benötigten Variablen
%%% p:v %%%
PlotVars.xscale = linspace(16.67,150,length(ds{Box_Cox,1}));
PlotVars.Xnew_1 = [PlotVars.mittel_a, PlotVars.mittel_f, linspace(16.67,150,length(ds{Box_Cox,1}))', PlotVars.p_1];
PlotVars.Xnew_2 = [PlotVars.mittel_a, PlotVars.mittel_f, linspace(16.67,150,length(ds{Box_Cox,1}))', PlotVars.p_2];
[PlotVars.ypred_1, PlotVars.yci_1] = predict(Models{ModelChoice,1}, PlotVars.Xnew_1, 'alpha', confi);
[PlotVars.ypred_2, PlotVars.yci_2] = predict(Models{ModelChoice,1}, PlotVars.Xnew_2, 'alpha', confi);
% Rücktrafo Box-Cox
if (Box_Cox == 1)
    PlotVars.ypred_1 = (lambda*(1/lambda + PlotVars.ypred_1)).^(1/lambda);
    PlotVars.ypred_2 = (lambda*(1/lambda + PlotVars.ypred_2)).^(1/lambda);
    PlotVars.yci_1 = (lambda*(1/lambda + PlotVars.yci_1)).^(1/lambda);
    PlotVars.yci_2 = (lambda*(1/lambda + PlotVars.yci_2)).^(1/lambda);
end
%%%
PlotVars.xLABEL = num2cell(50:50:150);
PlotVars.xLABEL = cellfun(@num2str, PlotVars.xLABEL, 'UniformOutput', 0);
PlotVars.xLABEL{1,length(PlotVars.xLABEL)-1} = 'mm/s';

subplot(4,4,15), plot(PlotVars.xscale, PlotVars.ypred_1*PlotVars.Umrechnung, 'Color', [1 0.5 0])
hold on
plot(PlotVars.xscale, PlotVars.ypred_2*PlotVars.Umrechnung, 'b')
plot(PlotVars.xscale, PlotVars.yci_1*PlotVars.Umrechnung, '--', 'Color', [1 0.5 0])
plot(PlotVars.xscale, PlotVars.yci_2*PlotVars.Umrechnung, '--b')
set(gca,'units','centimeters','Position',[8, 1, 3, 3], 'Fontsize', FontSize)
xlabel('v', 'Fontsize', FontSize+2)
set(gca, 'YTickLabel', {})
set(gca, 'XTick', 50:50:150)
set(gca, 'XTickLabel', PlotVars.xLABEL)
set(gca, 'YTick', 0:PlotVars.Step:PlotVars.IntyScale)
ylim([0 PlotVars.IntyScale])
xlim([16.67 150])
%%%
subplot(4,4,16), axis off

%% 3D Plots ohne Versuchspunkte
% Amplitude/Frequenz
PlotVars.xscale = linspace(0.2,0.75,length(ds{Box_Cox,1}))';
PlotVars.yscale = linspace(200,2000,length(ds{Box_Cox,1}))';

PlotVars.predictsZ = zeros(length(PlotVars.yscale), length(PlotVars.yscale));
for i = 1:length(PlotVars.xscale)
	PlotVars.predictsZ(:,i) = predict(Models{ModelChoice,1}, ...
        [repmat(PlotVars.xscale(i), [length(PlotVars.yscale) 1]), ...
        PlotVars.yscale, PlotVars.mittel_v, PlotVars.mittel_p]);
end

% Rücktrafo Box-Cox
if (Box_Cox == 1)
   PlotVars.predictsZ = arrayfun(box_back, PlotVars.predictsZ);
end

shades = linspace(0, 1, 50);
shades = repmat(shades, [3 1]);
shades = shades';

figure;
surf(PlotVars.xscale, PlotVars.yscale, PlotVars.predictsZ)
set(gca, 'YDir', 'reverse')
colormap(shades)
xlim([0.2, 0.75])
ylim([200, 2000])
xlabel('a')
ylabel('f')
zlabel('Zielgröße')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude/Vorschub
PlotVars.xscale = linspace(0.2,0.75,length(ds{Box_Cox,1}))';
PlotVars.yscale = linspace(16.67,150,length(ds{Box_Cox,1}))';

PlotVars.predictsZ = zeros(length(PlotVars.yscale), length(PlotVars.yscale));
for i = 1:length(PlotVars.xscale)
	PlotVars.predictsZ(:,i) = predict(Models{ModelChoice,1}, ... 
       [repmat(PlotVars.xscale(i), [length(PlotVars.yscale) 1]), ... 
       PlotVars.mittel_f, PlotVars.yscale, PlotVars.mittel_p]);
end

% Rücktrafo Box-Cox
if (Box_Cox == 1)
   PlotVars.predictsZ = arrayfun(box_back, PlotVars.predictsZ);
end

figure;
surf(PlotVars.xscale, PlotVars.yscale, PlotVars.predictsZ)
set(gca, 'YDir', 'reverse')
colormap(shades)
xlim([0.2, 0.75])
ylim([16.67, 150])
xlabel('a')
ylabel('v')
zlabel('Zielgröße')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Amplitude/Leistung
PlotVars.xscale = linspace(0.2,0.75,length(ds{Box_Cox,1}))';
PlotVars.yscale = linspace(1000,3000,length(ds{Box_Cox,1}))';

PlotVars.predictsZ = zeros(length(PlotVars.yscale), length(PlotVars.yscale));
for i = 1:length(PlotVars.xscale)
	PlotVars.predictsZ(:,i) = predict(Models{ModelChoice,1}, ...
        [repmat(PlotVars.xscale(i), [length(PlotVars.yscale) 1]), ...
        PlotVars.mittel_f, PlotVars.mittel_v, PlotVars.yscale]);
end

% Rücktrafo Box-Cox
if (Box_Cox == 1)
   PlotVars.predictsZ = arrayfun(box_back, PlotVars.predictsZ);
end

figure;
surf(PlotVars.xscale, PlotVars.yscale, PlotVars.predictsZ)
set(gca, 'YDir', 'reverse')
colormap(shades)
xlim([0.2, 0.75])
ylim([1000, 3000])
xlabel('a')
ylabel('p')
zlabel('Zielgröße')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Frequenz/Vorschub
PlotVars.xscale = linspace(200,2000,length(ds{Box_Cox,1}))';
PlotVars.yscale = linspace(16.67,150,length(ds{Box_Cox,1}))';

PlotVars.predictsZ = zeros(length(PlotVars.yscale), length(PlotVars.yscale));
for i = 1:length(PlotVars.xscale)
	PlotVars.predictsZ(:,i) = predict(Models{ModelChoice,1}, ...
        [PlotVars.mittel_a, repmat(PlotVars.xscale(i), ...
        [length(PlotVars.yscale) 1]), PlotVars.yscale, PlotVars.mittel_p]);
end

% Rücktrafo Box-Cox
if (Box_Cox == 1)
   PlotVars.predictsZ = arrayfun(box_back, PlotVars.predictsZ);
end

figure;
surf(PlotVars.xscale, PlotVars.yscale, PlotVars.predictsZ)
set(gca, 'YDir', 'reverse')
colormap(shades)
xlim([200, 2000])
ylim([16.67, 150])
xlabel('f')
ylabel('v')
zlabel('Zielgröße')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Frequenz/Leistung
PlotVars.xscale = linspace(200,2000,length(ds{Box_Cox,1}))';
PlotVars.yscale = linspace(1000,3000,length(ds{Box_Cox,1}))';

PlotVars.predictsZ = zeros(length(PlotVars.yscale), length(PlotVars.yscale));
for i = 1:length(PlotVars.xscale)
	PlotVars.predictsZ(:,i) = predict(Models{ModelChoice,1}, ...
        [PlotVars.mittel_a, repmat(PlotVars.xscale(i), ...
        [length(PlotVars.yscale) 1]), PlotVars.mittel_v, PlotVars.yscale]);
end

% Rücktrafo Box-Cox
if (Box_Cox == 1)
   PlotVars.predictsZ = arrayfun(box_back, PlotVars.predictsZ);
end

figure;
surf(PlotVars.xscale, PlotVars.yscale, PlotVars.predictsZ)
set(gca, 'YDir', 'reverse')
colormap(shades)
xlim([200, 2000])
ylim([1000, 3000])
xlabel('f')
ylabel('p')
zlabel('Zielgröße')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vorschub/Leistung
PlotVars.xscale = linspace(16.67,150,length(ds{Box_Cox,1}))';
PlotVars.yscale = linspace(1000,3000,length(ds{Box_Cox,1}))';

PlotVars.predictsZ = zeros(length(PlotVars.yscale), length(PlotVars.yscale));
for i = 1:length(PlotVars.xscale)
	PlotVars.predictsZ(:,i) = predict(Models{ModelChoice,1}, ...
        [PlotVars.mittel_a, PlotVars.mittel_f, repmat(PlotVars.xscale(i), ...
        [length(PlotVars.yscale) 1]), PlotVars.yscale]);
end

% Rücktrafo Box-Cox
if (Box_Cox == 1)
   PlotVars.predictsZ = arrayfun(box_back, PlotVars.predictsZ);
end

figure;
surf(PlotVars.xscale, PlotVars.yscale, PlotVars.predictsZ)
set(gca, 'YDir', 'reverse')
colormap(shades)
xlim([16.67, 150])
ylim([1000, 3000])
xlabel('v')
ylabel('p')
zlabel('Zielgröße')