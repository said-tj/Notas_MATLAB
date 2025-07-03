function PolyVIEW_Replica
% PolyVIEW_Replica  Programa para replicar PolyVIEW en MATLAB
%   Carga el .BIN y .LOG, aplica escala a unidades físicas y ofrece
%   una GUI con menú desplegable de canal y slider de duración.

    %% — 1. Definición de rutas y parámetros
    binFile = 'data/ratpiloto20092024b.BIN';
    logFile = 'data/ratpiloto20092024b.LOG';
    fs      = 300;                          % frecuencia de muestreo [Hz]

    %% — 2. Carga de datos
    [data, meta] = loadPolyViewBIN(binFile);
    events       = loadPolyViewLog(logFile);

    t = (0:size(data,1)-1)/fs;              % eje temporal en segundos

    %% — 3. Creación de la figura y controles
    hFig   = figure('Name','PolyVIEW Replica', ...
                    'NumberTitle','off', ...
                    'Position',[100 100 900 600]);

    hAx    = axes('Parent',hFig, ...
                  'Units','normalized', ...
                  'Position',[0.08 0.30 0.88 0.65]);

    % Dropdown para seleccionar canal
    hPopup = uicontrol('Style','popupmenu', ...
                       'Parent',hFig, ...
                       'Units','normalized', ...
                       'Position',[0.08 0.18 0.30 0.05], ...
                       'String',{meta.name}, ...
                       'Callback',@updatePlot);

    % Slider para la duración de la ventana (en segundos)
    hSlider = uicontrol('Style','slider', ...
                        'Parent',hFig, ...
                        'Units','normalized', ...
                        'Position',[0.45 0.18 0.45 0.05], ...
                        'Min',1,'Max',120,'Value',10, ...
                        'SliderStep',[1/119 10/119], ...
                        'Callback',@updatePlot);

    % Texto que muestra la duración seleccionada
    hText = uicontrol('Style','text', ...
                      'Parent',hFig, ...
                      'Units','normalized', ...
                      'Position',[0.08 0.10 0.30 0.05], ...
                      'String','Duración [s]: 10');

    %% — 4. Primer dibujo
    updatePlot()

    %% ──────────────────────────────────────────────────────────────────────
    function updatePlot(~,~)
    % Callback que redibuja la señal según canal y duración
        chIdx = get(hPopup,'Value');                % índice de canal
        dur   = round(get(hSlider,'Value'));        % segundos a mostrar
        set(hText,'String',sprintf('Duración [s]: %d',dur))

        N     = dur * fs;                           % muestras a mostrar
        i0    = 1;                                  % siempre inicio en el sample 1
        i1    = min(size(data,1), N);

        % Graficar señal escalada
        plot(hAx, t(i0:i1), data(i0:i1, chIdx), 'LineWidth',1.2)
        hold(hAx,'on')

        % Superponer líneas de evento y etiqueta
        selE = events.sample <= i1;
        for k = find(selE)'
            xx = events.sample(k)/fs;
            yy = interp1(t(i0:i1), data(i0:i1,chIdx), xx);
            line(hAx, [xx xx], hAx.YLim, 'Color','r','LineStyle','--')
            text(hAx, xx, yy, events.label{k}, ...
                 'VerticalAlignment','bottom','Color','r','FontSize',8)
        end

        hold(hAx,'off')
        xlabel(hAx,'Tiempo [s]')
        ylabel(hAx, meta(chIdx).unit)
        title(hAx, sprintf('%s  (%s)', meta(chIdx).name, meta(chIdx).unit))
        grid(hAx,'on')
    end
end

%% ────────────────────────────────────────────────────────────────────────────
function [data, meta] = loadPolyViewBIN(filename)
% loadPolyViewBIN  Lee .BIN de PolyVIEW y devuelve señal + metadata
    fid     = fopen(filename,'rb');
    hdrSize = fread(fid,1,'uint32','ieee-be');    % tamaño cabecera
    nCh     = fread(fid,1,'uint32','ieee-be');    % nº de canales

    % Leemos la cabecera completa (aquí podrías parsear, pero lo hacemos
    % manualmente según CTBC)
    fseek(fid, hdrSize, 'bof');
    fread(fid, 4, 'uint8');                        % padding 4 bytes

    % Leemos todas las muestras como int16 little-endian
    raw = fread(fid, inf, 'int16=>double', 'ieee-le');
    fclose(fid);

    % Reorganizamos: cada fila = instante, cada columna = canal
    data = reshape(raw, nCh, []).';

    % Metadatos manuales (ajusta según tu calibración real)
    meta(1).name   = 'Ch1'; meta(1).unit   = 'V';    meta(1).scale = 3.0519e-4; meta(1).offset = 0; meta(1).gain = 1;
    meta(2).name   = 'Ch2'; meta(2).unit   = 'V';    meta(2).scale = 3.0519e-4; meta(2).offset = 0; meta(2).gain = 1;
    meta(3).name   = 'Ch3'; meta(3).unit   = 'mmHg'; meta(3).scale = 6.1037e-2; meta(3).offset = 0; meta(3).gain = 1;

    % Aplicamos escala → unidades físicas
    for k = 1:nCh
        data(:,k) = data(:,k) * meta(k).scale * meta(k).gain + meta(k).offset;
    end
end

%% ────────────────────────────────────────────────────────────────────────────
function events = loadPolyViewLog(filename)
% loadPolyViewLog  Lee .LOG de PolyVIEW y devuelve tabla de eventos
    fid = fopen(filename,'r');
    C   = textscan(fid, '%s %f %f', 'CommentStyle','*');
    fclose(fid);
    events = table(C{1}, C{2}, C{3}, 'VariableNames', {'label','sample','code'});
end
