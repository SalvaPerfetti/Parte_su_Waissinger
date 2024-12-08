close all
clear all
clc

%% Test Case 1

U_Inf_Mag = 1;
beta = 30;
U_Inf = [cosd(beta) sind(beta) 0] .* U_Inf_Mag;
rho = 1.225;

config.NCorpi = 1; % Numero di "corpi"

config.RootChord = [1]; %: Corda alla radice dell'ala
config.DihedralAngle = [30]; % [°] %Angolo di diedro dell'ala 
config.SweepAngle = [1]; % [°] %Angolo di freccia 
config.TaperRatio = [1]; % rapporto di rastremazione: Corda radice/ Corda all'estremità
config.AspectRatio = [5]; % b^2/S
config.Span = 2; % B:apertura alare
config.LEPosition_X = [0]; % Coordinata x della posizione del bordo d'attacco (Leading Edge, LE) della radice dell'ala.
config.LEPosition_Y = [0]; % Coordinata y della posizione del bordo d'attacco della radice.
config.LEPosition_Z = [0]; % Coordinata z della posizione del bordo d'attacco della radice.

config.RotationAngle_X = [0]; % Angolo di rotazione attorno all'asse x, permette di inclinare l'ala nel piano yz
config.RotationAngle_Y = [4]; % Angolo di rotazione attorno all'asse y
config.RotationAngle_Z = [0]; % Angolo di rotazione attorno all'asse z 

% Discretization options
config.SemiSpanwiseDiscr = [30]; % pannelli lungo la semi-apertura (b/2)
config.ChordwiseDiscr = [20]; % pannelli lungo la corda


%% Preliminary computations

% Computing the span (Calcola la semi-apertura (b/2) dell'ala)
config.SemiSpan = config.Span./2; 
% Computing the surface (Calcola la superficie alare totale (S))
config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2);
% Calcola la superficie alare proiettata sul piano orizzontale, tenendo conto del diedro
config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
% Computing the Tip chord (Calcola la corda all' estremità alare)
config.TipChord = config.RootChord .* config.TaperRatio;

% Compute MAC (Calcola la corda media aerodinamica)
config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));

%% Create the geometry structure

ControlPoints = cell(config.NCorpi, 1); % punti di controllo dove vengono calcolati gli effetti aerodinamici dei vortici.
InducedPoints = cell(config.NCorpi, 1); % punti intermedi lungo il vortice usati per calcolare l'induzione
Normals = cell(config.NCorpi, 1); % vettori normali ai pannelli, utilizzati per calcolare le influenze
InfiniteVortices = cell(config.NCorpi, 1); % coordinate dei vortici semi-infiniti
Vortices = cell(config.NCorpi, 1); % coordinate dei vortici finiti (da radice a estremità)
internalMesh = cell(config.NCorpi, 1); % mesh interna dell'ala
WingExtremes = cell(config.NCorpi, 1); % coordinate degli estremi dell'ala (bordo d'attacco e bordo d'uscita sia alla radice che alla punta).


for iCorpo = 1:config.NCorpi

    [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure(config, iCorpo);

end
    

%% Matrices initialization

NPanelsTot = 2* config.SemiSpanwiseDiscr * config.ChordwiseDiscr'; % Calcola il numero totale di pannelli sulla superficie dell'ala
matriceA = zeros(NPanelsTot, NPanelsTot); % Matrice quadrata di dimensioni NPanelsTot×NPanelsTot
TermineNoto = zeros(NPanelsTot, 1); % Vettore colonna di dimensione NPanelsTot

%% Construction of the matrix

rowIndex = 0;
for iCorpo = 1:config.NCorpi 
    
    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo) % Indice per i pannelli lungo la corda 
        % Cycle on all of its spanwise panels (lungo l'apertura) 
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
   
            columnIndex = 0;
            
            ControlPointHere = ControlPoints{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords; % Coordinate del punto di controllo del pannello attuale (centro del pannello).
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords; %  Vettore normale al pannello corrente.
            
            
            for jCorpo = 1:config.NCorpi % Corpo che induce l’influenza.
                
                % Cycle on all of its chordwise panels
                for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo) %Indici del pannello che influenza il pannello corrente
                    % Cycle on all of its spanwise panels
                    for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo) %Indici del pannello che influenza il pannello corrente
                        
                        % Update column index
                        columnIndex = columnIndex + 1;
                        
                        % Compute the influence induced by first
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                        U = vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);
% Calcola il contributo del primo vortice semi-infinito (radice del pannello che si estende verso l'infinito).
                        
                        % Compute the influence induced by finite vortex
                        Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                        Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);
% Calcola il contributo del vortice finito tra i due estremi del pannello influente (radice e punta).
                        
                        % Compute the influence induced by second
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);
% Calcola il contributo del secondo vortice semi-infinito (punta del pannello che si estende verso l'infinito).
                        
                        matriceA(rowIndex, columnIndex) = dot(U, NormalHere); %Prodotto scalare tra l’influenza U e la normale del pannello corrente. Rappresenta l'effetto aerodinamico del pannello j sul pannello i.
                        
                    end
                end

            end
            
        
            
        end
    end
end

%% Costruzione del termine noto
rowIndex = 0;
for iCorpo = 1:config.NCorpi
    
    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
  
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            
            TermineNoto(rowIndex) = -dot(U_Inf, NormalHere);
            
        end
    end
end

%% Solve the linear system

Solution = linsolve(matriceA, TermineNoto);

Gamma = cell(config.NCorpi, 1);

rowIndex = 0;
for iCorpo = 1:config.NCorpi
    
    Gamma{iCorpo} = zeros( config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 ); % definisco Gamma
    
     % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
            
            Gamma{iCorpo}(ChordPanel_i, SpanPanel_i) = Solution(rowIndex);
        end
        
    end
    
end

%% Visualization
% % PUNTO 1: Visualizzazione della circolazione
% figure;
% for iCorpo = 1:config.NCorpi
%     surf(1:config.SemiSpanwiseDiscr(iCorpo)*2, 1:config.ChordwiseDiscr(iCorpo), Gamma{iCorpo});
%     xlabel('Pannelli in apertura');
%     ylabel('Pannelli in corda');
%     zlabel('\Gamma (Circolazione)');
%     title('Distribuzione della circolazione sui pannelli');
%     colorbar;
%     view(3);
% end
% 

%% Visualization of Circulation Distribution
figure;
for iCorpo = 1:config.NCorpi
    [X, Y] = meshgrid(1:config.SemiSpanwiseDiscr(iCorpo)*2, 1:config.ChordwiseDiscr(iCorpo));
    surf(X, Y, Gamma{iCorpo});
    xlabel('Spanwise Panels');
    ylabel('Chordwise Panels');
    zlabel('\Gamma [Circulation]');
    title('Circulation Distribution');
    colormap jet;
    colorbar;
    hold on;
end
hold off;

%% Compute the 2D and 3D Lift
% % PUNTO 2: Calcolo della portanza 
% Lift_2D = cell(config.NCorpi, 1);
% Lift_3D = 0;
% 
% for iCorpo = 1:config.NCorpi
%     SemiSpanDiscr = config.SemiSpanwiseDiscr(iCorpo);
%     Delta_s = config.Span / (2 * SemiSpanDiscr); % Lunghezza del pannello in apertura
%     Lift_2D{iCorpo} = rho * Gamma{iCorpo} * Delta_s; % Lift 2D
%     Lift_3D = Lift_3D + sum(Lift_2D{iCorpo}, 'all'); % Lift totale 3D
% end
% 
% disp(['Lift totale (3D): ', num2str(Lift_3D), ' N']);

%% Lift Computation
Lift2D = cell(config.NCorpi, 1); % Portanza distribuita lungo l'apertura alare
TotalLift = 0;

for iCorpo = 1:config.NCorpi
    Lift2D{iCorpo} = zeros(1, config.SemiSpanwiseDiscr(iCorpo)*2);
    dSpan = config.Span / (config.SemiSpanwiseDiscr(iCorpo)*2); % Discretizzazione in apertura
    
    for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
        Lift2D{iCorpo}(SpanPanel_i) = rho * U_Inf_Mag * sum(Gamma{iCorpo}(:, SpanPanel_i)) ; % * dSpan
    end
    
    TotalLift = TotalLift + sum(Lift2D{iCorpo});
end

fprintf('Total Lift: %.3f N\n', TotalLift);

% Visualizzazione della portanza distribuita
figure;
for iCorpo = 1:config.NCorpi
    plot(1:config.SemiSpanwiseDiscr(iCorpo)*2, Lift2D{iCorpo}, 'LineWidth', 1.5);
    hold on;
end
xlabel('Spanwise Panels');
ylabel('Lift [N]');
title('Distributed Lift Along the Span');
legend('Corpo 1');
grid on;
hold off;

%% Compute 2D and 3D induced drag
% % PUNTO 3: Calcolo della resistenza indotta
% e = 0.9; % Fattore di efficienza
% AR = config.AspectRatio; % Aspect Ratio
% Induced_Drag_3D = (Lift_3D^2) / (pi * e * AR);
% 
% disp(['Resistenza indotta (3D): ', num2str(Induced_Drag_3D), ' N']);
%% Induced Drag Computation
InducedDrag = 0;

for iCorpo = 1:config.NCorpi
    for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
        vind = 0;
        
        % Calcola l'induzione dal wake (vortici semi-infiniti)
        for ChordPanel_j = 1:config.ChordwiseDiscr(iCorpo)
            Extreme_1 = Vortices{iCorpo}{ChordPanel_j, SpanPanel_i}.Root;
            Extreme_2 = Vortices{iCorpo}{ChordPanel_j, SpanPanel_i}.Tip;
            ControlPointHere = ControlPoints{iCorpo}{ChordPanel_j, SpanPanel_i}.Coords;
            vind = vind + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);
        end
        
        % Calcolo angolo di incidenza indotto
        alphaInd = atan2(norm(vind), U_Inf_Mag);
        
        % Calcolo della resistenza indotta
        LiftHere = rho * U_Inf_Mag * sum(Gamma{iCorpo}(:, SpanPanel_i)) * dSpan;
        InducedDrag = InducedDrag + LiftHere * sin(alphaInd);
    end
end

fprintf('Induced Drag: %.3f N\n', InducedDrag);

%% Grafico della distribuzione di portanza sull'ala
%% Lift Distribution Along the Span
figure;
for iCorpo = 1:config.NCorpi
    % Calcolo delle posizioni spanwise (centri dei pannelli)
    y = linspace(-config.Span/2, config.Span/2, 2*config.SemiSpanwiseDiscr(iCorpo));
    LiftDistribution = sum(Gamma{iCorpo}, 1) * rho * U_Inf_Mag; % Distribuzione della portanza
    
    plot(y, LiftDistribution, 'LineWidth', 1.5);
    hold on;
end

xlabel('Spanwise Position [m]');
ylabel('Lift Distribution [N/m]');
title('Lift Distribution Along the Wing Span');
legend('Corpo 1');
grid on;
hold off;

