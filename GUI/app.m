%--- Copyright notice ---%
% Copyright (C) 2026 The University of Manchester
% Written by David Mostaza Prieto,  Nicholas H. Crisp, Luciana Sinpetru, 
% Sabrina Livadiotti and Joseph Tucker
%
% This file is part of the ADBSat toolkit.
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%------------ BEGIN CODE ----------%
classdef app < matlab.apps.AppBase
    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        MainGridLayout                  matlab.ui.container.GridLayout
        LeftPanel                       matlab.ui.container.Panel
        TabGroup                        matlab.ui.container.TabGroup
        EnvironmentTab                  matlab.ui.container.Tab
        APMagneticIndexEditField_2      matlab.ui.control.EditField
        EditFieldLabel                  matlab.ui.control.Label
        AnomalousOxygenCheckBox         matlab.ui.control.CheckBox
        F107dailyEditField              matlab.ui.control.NumericEditField
        F107dailyEditFieldLabel         matlab.ui.control.Label
        F107averageEditField            matlab.ui.control.NumericEditField
        F107averageEditFieldLabel       matlab.ui.control.Label
        TimeofdaysecondsEditField       matlab.ui.control.NumericEditField
        TimeofdaysecondsEditFieldLabel  matlab.ui.control.Label
        DayofYearDatePicker             matlab.ui.control.DatePicker
        DayofYearDatePickerLabel        matlab.ui.control.Label
        LongitudedegEditField           matlab.ui.control.NumericEditField
        LongitudedegEditFieldLabel      matlab.ui.control.Label
        InclinationdegEditField         matlab.ui.control.NumericEditField
        InclinationdegEditFieldLabel    matlab.ui.control.Label
        AltitudekmEditField             matlab.ui.control.NumericEditField
        AltitudekmEditFieldLabel        matlab.ui.control.Label
        SimulationTab                   matlab.ui.container.Tab
        PlotDropDown                    matlab.ui.control.DropDown
        PlotDropDownLabel               matlab.ui.control.Label
        RunButton                       matlab.ui.control.Button
        CancelButton                    matlab.ui.control.Button
        ClearPlotsButton                matlab.ui.control.Button
        delCheckBox                     matlab.ui.control.CheckBox
        VerboseCheckBox                 matlab.ui.control.CheckBox
        DiffuseReflecitvityEditField    matlab.ui.control.NumericEditField
        EditField_2Label_2              matlab.ui.control.Label
        SpecularReflectivityEditField   matlab.ui.control.NumericEditField
        EditField_2Label                matlab.ui.control.Label
        SolarRadiationPressureCheckBox  matlab.ui.control.CheckBox
        WallTermperatureEditField       matlab.ui.control.EditField
        WallTermperatureKelvinLabel     matlab.ui.control.Label
        EditField                       matlab.ui.control.EditField
        AccomodationCoefficientsLabel   matlab.ui.control.Label
        AngleofSideslipEditField        matlab.ui.control.NumericEditField
        AngleofSideslipEditFieldLabel   matlab.ui.control.Label
        AngleofAttackEditField          matlab.ui.control.NumericEditField
        AngleofAttackEditFieldLabel     matlab.ui.control.Label
        HyperthermalCheckBox            matlab.ui.control.CheckBox
        RaysEditField                   matlab.ui.control.NumericEditField
        RaysEditFieldLabel              matlab.ui.control.Label
        ShadowCheckBox                  matlab.ui.control.CheckBox
        GSIModelDropDown                matlab.ui.control.DropDown
        GSIModelDropDownLabel           matlab.ui.control.Label
        ModelTab                        matlab.ui.container.Tab
        PitchEditField                  matlab.ui.control.NumericEditField
        PitchEditFieldLabel             matlab.ui.control.Label
        YawEditField                    matlab.ui.control.NumericEditField
        YawEditFieldLabel               matlab.ui.control.Label
        RollEditField                   matlab.ui.control.NumericEditField
        RollLabel                       matlab.ui.control.Label
        RotationLabel                   matlab.ui.control.Label
        TranslationLabel                matlab.ui.control.Label
        inZEditField                    matlab.ui.control.NumericEditField
        inZLabel                        matlab.ui.control.Label
        inYEditField                    matlab.ui.control.NumericEditField
        inYLabel                        matlab.ui.control.Label
        inXEditField                    matlab.ui.control.NumericEditField
        EditField_3Label_3              matlab.ui.control.Label
        AutoCenteringCheckBox           matlab.ui.control.CheckBox
        ElementSizemEditField           matlab.ui.control.EditField
        EditField_3Label_2              matlab.ui.control.Label
        NumberofSubdivisionsEditField   matlab.ui.control.NumericEditField
        EditField_3Label                matlab.ui.control.Label
        
        % Model Selection UI
        BrowseButton                    matlab.ui.control.Button
        SelectedPathEditField           matlab.ui.control.EditField
        SelectedPathLabel               matlab.ui.control.Label
        ProcessGeometryButton           matlab.ui.control.Button
        
        % Preset save/load
        SavePresetButton                matlab.ui.control.Button
        LoadPresetButton                matlab.ui.control.Button
        
        % Export controls
        ExportPlotButton                matlab.ui.control.Button
        ExportResultsButton             matlab.ui.control.Button
        ExportAllPlotsButton            matlab.ui.control.Button

        % Plot display controls
        LockColorScaleCheckBox          matlab.ui.control.CheckBox
        
        % Right Panel Structure
        RightPanel                      matlab.ui.container.Panel
        RightTabGroup                   matlab.ui.container.TabGroup
        PlotTab1                        matlab.ui.container.Tab
        PlotTab2                        matlab.ui.container.Tab
        PlotTab4                        matlab.ui.container.Tab
        PlotTab5                        matlab.ui.container.Tab
        PlotTab6                        matlab.ui.container.Tab
        PlotTab7                        matlab.ui.container.Tab
        UIAxes1                         matlab.ui.control.UIAxes
        UIAxes2                         matlab.ui.control.UIAxes
        UIAxesMeshFinal                 matlab.ui.control.UIAxes
        MeshViewSwitch                  matlab.ui.control.Switch
        UIAxesQuality                   matlab.ui.control.UIAxes
        UIAxesNormals                   matlab.ui.control.UIAxes
        UIAxesMaterial                  matlab.ui.control.UIAxes
        ResultsTable                    matlab.ui.control.Table
        
        % Bottom Log Console Structural Elements
        BottomPanel                     matlab.ui.container.Panel
        ConsoleTextArea                 matlab.ui.control.TextArea
        
        % Progress Bar specific to bottom left - does not work yet
        ProgressBar                     matlab.ui.control.LinearGauge
        
        % Target internal paths extracted globally from Browse workflow
        UserSelectedFullFile            char = ''
    end
    
    properties (Access = private)
        onePanelWidth = 640;
        CancelRequested logical = false; % flag to handle loop escapes
        LastFileOut char = '';          % Path to most recent aero results .mat
        LastModOut char = '';           % Path to most recent mesh .mat
        LastSummaryTable table = table.empty(0,2); % Cached results-summary table for CSV export
        LockedCLim double = []; % Captured colour range while Lock Color Scale is on
    end
    
    methods (Access = public)
        % Print tracking string arrays directly to terminal screen
        function logMessage(app, msgText)
            currentText = app.ConsoleTextArea.Value;
            if isempty(currentText) || (numel(currentText) == 1 && isempty(currentText{1}))
                app.ConsoleTextArea.Value = {msgText};
            else
                app.ConsoleTextArea.Value = [currentText; {msgText}];
            end
            scroll(app.ConsoleTextArea, 'bottom');
            drawnow;
        end
        
        function clearLog(app)
            app.ConsoleTextArea.Value = '';
        end
        
        function clearAllPlots(app)
            axesList = [app.UIAxes1, app.UIAxes2, app.UIAxesMeshFinal, ...
                        app.UIAxesQuality, app.UIAxesNormals, app.UIAxesMaterial];
            for k = 1:numel(axesList)
                cla(axesList(k));
                legend(axesList(k), 'off');
                colorbar(axesList(k), 'off');
                title(axesList(k), '');
                subtitle(axesList(k), '');
            end
            app.logMessage('Cleared all plot axes.');
        end
        
        % Draws the 3D surface-response plot (UIAxes1) for a given cached
        % results file / mesh, coloured by whichever field is passed in.
        % Factored out so it can be called both right after a solve and
        % again later purely to switch the displayed field.
        function renderSurfacePlot(app, fileOut, modOut, aoa_deg, aos_deg, plotField)
            app.clearAerodynamicPlot();

            plot_surfq_app(fileOut, modOut, aoa_deg, aos_deg, plotField, app.UIAxes1);

            cb = colorbar(app.UIAxes1);
            cb.Label.String = upper(plotField);

            title(app.UIAxes1, sprintf('Surface Distribution [%s]', upper(plotField)), 'FontSize', 12);
            subtitle(app.UIAxes1, sprintf('AoA (\\alpha): %.1f°  |  AoS (\\beta): %.1f°', aoa_deg, aos_deg), 'FontSize', 10);

            patches = findobj(app.UIAxes1, 'Type', 'Patch');
            for k = 1:numel(patches)
                set(patches(k), 'EdgeColor', 'none');
            end

            app.applyColorScaleLock();
        end

        % If "Lock Color Scale" is on, re-applies whatever colour range was
        % captured when it was switched on, so later redraws (a different
        % field, a different angle, or a fresh run) stay visually
        % comparable instead of each one auto-scaling to its own data.
        function applyColorScaleLock(app)
            if app.LockColorScaleCheckBox.Value
                if isempty(app.LockedCLim)
                    app.LockedCLim = app.UIAxes1.CLim;
                else
                    app.UIAxes1.CLim = app.LockedCLim;
                end
            end
        end

        % Captures the currently displayed colour range as a fixed scale
        % when turned on, or clears it and resumes auto-scaling when
        % turned off.
        function LockColorScaleCheckBoxValueChanged(app, event)
            if app.LockColorScaleCheckBox.Value
                app.LockedCLim = app.UIAxes1.CLim;
                app.logMessage(sprintf('Color scale locked to [%.3g, %.3g].', app.LockedCLim(1), app.LockedCLim(2)));
            else
                app.LockedCLim = [];
                app.logMessage('Color scale unlocked (auto-scaling resumed).');
            end
        end

        % Shared logic for re-rendering the cached surface-response plot
        % with whichever field is currently selected in PlotDropDown,
        % without re-running the (potentially expensive) solver pipeline.
        function refreshSurfacePlot(app, notifyIfMissing)
            if nargin < 2
                notifyIfMissing = false;
            end
            if isempty(app.LastFileOut) || ~exist(app.LastFileOut, 'file') || ...
                    isempty(app.LastModOut) || ~exist(app.LastModOut, 'file')
                if notifyIfMissing
                    app.logMessage('No aerodynamic results are loaded yet - run a simulation first.');
                end
                return;
            end

            aoa_deg = app.AngleofAttackEditField.Value;
            aos_deg = app.AngleofSideslipEditField.Value;

            try
                app.renderSurfacePlot(app.LastFileOut, app.LastModOut, aoa_deg(1), aos_deg(1), app.PlotDropDown.Value);
                app.logMessage(sprintf('Plot updated: %s', app.PlotDropDown.Value));
            catch ME
                app.logMessage(sprintf('PLOT FAILED: %s', ME.message));
            end
        end

        % Lets the user switch which result field is shown on the 3D
        % surface without re-running the simulation, as long as results
        % from a previous run are still cached.
        function PlotDropDownValueChanged(app, event)
            app.refreshSurfacePlot(false);
        end

        % Clears only the aerodynamic surface-response plot (UIAxes1).

        % Used by the run sim button as all other plots are now produced by
        % the import mesh control
        function clearAerodynamicPlot(app)
            cla(app.UIAxes1);
            legend(app.UIAxes1, 'off');
            colorbar(app.UIAxes1, 'off');
            title(app.UIAxes1, '');
            subtitle(app.UIAxes1, '');
        end
        
        % Shows/hides one of the overlaid mesh axes this is for the mesh
        % toggle plot
        function setMeshAxesVisibility(~, ax, tf)
            if tf
                ax.Visible = 'on';
                set(ax.Children, 'Visible', 'on');
            else
                ax.Visible = 'off';
                set(ax.Children, 'Visible', 'off');
            end
        end
        
        % Applies the Mesh view switch's current selection to both
        % overlaid mesh axes. Called on toggle, and again after every
        % import so newly (re)plotted content respects the current choice.
        function applyMeshViewToggle(app)
            showFinal = strcmp(app.MeshViewSwitch.Value, 'Final');
            app.setMeshAxesVisibility(app.UIAxes2, ~showFinal);
            app.setMeshAxesVisibility(app.UIAxesMeshFinal, showFinal);
        end
    end
    
    methods (Access = private)
        % open file window when browse button pushed
        function BrowseButtonPushed(app, event)
            [file, path] = uigetfile({ ...
                '*.obj;*.stp;*.stl', 'Supported Geometry Formats (*.obj, *.stp, *.stl)'; ...
                '*.obj', 'Wavefront Object Files (*.obj)'; ...
                '*.stp', 'STEP Exchange Files (*.stp)'; ...
                '*.stl', 'Stereolithography Files (*.stl)'}, ...
                'Select Satellite Geometry Mesh Model');
            
            if isequal(file, 0) || isequal(path, 0)
                app.logMessage('File selection cancelled by user.');
            else
                app.UserSelectedFullFile = fullfile(path, file);
                app.SelectedPathEditField.Value = app.UserSelectedFullFile;
                app.logMessage(sprintf('Geometry linked: %s', file));
            end
        end

        % Manual reset of every plot tab back to a blank state - likely
        % redundant now but previously useful
        function ClearPlotsButtonPushed(app, event)
            app.clearAllPlots();
        end
        
        % Toggle which of the overlaid Original/Final mesh axes is shown
        function MeshViewSwitchValueChanged(app, event)
            app.applyMeshViewToggle();
        end
        
        % Save every configuration to a JSON preset file
        function SavePresetButtonPushed(app, event)
            cfg = struct();
            cfg.presetVersion = 1;
            cfg.modelFile = app.UserSelectedFullFile;
            
            % Model Geometry tab
            cfg.subdivisions = app.NumberofSubdivisionsEditField.Value;
            cfg.elementSize = app.ElementSizemEditField.Value;
            cfg.autoCentering = app.AutoCenteringCheckBox.Value;
            cfg.translation = [app.inXEditField.Value, app.inYEditField.Value, app.inZEditField.Value];
            cfg.rotation = [app.PitchEditField.Value, app.RollEditField.Value, app.YawEditField.Value];
            
            % Environment tab
            cfg.altitude = app.AltitudekmEditField.Value;
            cfg.inclination = app.InclinationdegEditField.Value;
            cfg.longitude = app.LongitudedegEditField.Value;
            dpVal = app.DayofYearDatePicker.Value;
            if ~isempty(dpVal)
                cfg.dayDate = char(string(dpVal, 'yyyy-MM-dd'));
            else
                cfg.dayDate = '';
            end
            cfg.timeOfDay = app.TimeofdaysecondsEditField.Value;
            cfg.f107average = app.F107averageEditField.Value;
            cfg.f107daily = app.F107dailyEditField.Value;
            cfg.apIndex = app.APMagneticIndexEditField_2.Value;
            cfg.anomalousOxygen = app.AnomalousOxygenCheckBox.Value;
            
            % Simulation tab
            cfg.aoa = app.AngleofAttackEditField.Value;
            cfg.aos = app.AngleofSideslipEditField.Value;
            cfg.gsiModel = app.GSIModelDropDown.Value;
            cfg.accommodation = app.EditField.Value;
            cfg.wallTemp = app.WallTermperatureEditField.Value;
            cfg.shadow = app.ShadowCheckBox.Value;
            cfg.hyperthermal = app.HyperthermalCheckBox.Value;
            cfg.rays = app.RaysEditField.Value;
            cfg.solarRadiation = app.SolarRadiationPressureCheckBox.Value;
            cfg.specularReflectivity = app.SpecularReflectivityEditField.Value;
            cfg.diffuseReflectivity = app.DiffuseReflecitvityEditField.Value;
            cfg.verbose = app.VerboseCheckBox.Value;
            cfg.deleteTemp = app.delCheckBox.Value;
            cfg.plotField = app.PlotDropDown.Value;
            cfg.lockColorScale = app.LockColorScaleCheckBox.Value;
            
            [file, path] = uiputfile('*.json', 'Save Parameter Preset As', 'adbsat_preset.json');
            if isequal(file, 0)
                return;
            end
            
            try
                jsonStr = jsonencode(cfg);
                fid = fopen(fullfile(path, file), 'w');
                if fid == -1
                    error('Could not open file for writing.');
                end
                fwrite(fid, jsonStr, 'char');
                fclose(fid);
                app.logMessage(sprintf('Preset saved: %s', file));
            catch ME
                app.logMessage(sprintf('PRESET SAVE FAILED: %s', ME.message));
            end
        end
        
        % Load a JSON preset file and reassign inputs
        function LoadPresetButtonPushed(app, event)
            [file, path] = uigetfile('*.json', 'Load Parameter Preset');
            if isequal(file, 0)
                return;
            end
            
            try
                raw = fileread(fullfile(path, file));
                cfg = jsondecode(raw);
            catch ME
                app.logMessage(sprintf('PRESET LOAD FAILED: %s', ME.message));
                return;
            end
            
            if isfield(cfg, 'modelFile') && ~isempty(cfg.modelFile)
                if exist(cfg.modelFile, 'file')
                    app.UserSelectedFullFile = cfg.modelFile;
                    app.SelectedPathEditField.Value = cfg.modelFile;
                else
                    app.logMessage(sprintf('Preset model file not found on disk (please re-Browse): %s', cfg.modelFile));
                end
            end
            
            if isfield(cfg,'subdivisions'),   app.NumberofSubdivisionsEditField.Value = cfg.subdivisions;   end
            if isfield(cfg,'elementSize'),     app.ElementSizemEditField.Value = cfg.elementSize;             end
            if isfield(cfg,'autoCentering'),   app.AutoCenteringCheckBox.Value = logical(cfg.autoCentering);  end
            if isfield(cfg,'translation') && numel(cfg.translation) == 3
                app.inXEditField.Value = cfg.translation(1);
                app.inYEditField.Value = cfg.translation(2);
                app.inZEditField.Value = cfg.translation(3);
            end
            if isfield(cfg,'rotation') && numel(cfg.rotation) == 3
                app.PitchEditField.Value = cfg.rotation(1);
                app.RollEditField.Value  = cfg.rotation(2);
                app.YawEditField.Value   = cfg.rotation(3);
            end
            
            if isfield(cfg,'altitude'),     app.AltitudekmEditField.Value = cfg.altitude;         end
            if isfield(cfg,'inclination'),  app.InclinationdegEditField.Value = cfg.inclination;  end
            if isfield(cfg,'longitude'),    app.LongitudedegEditField.Value = cfg.longitude;      end
            if isfield(cfg,'dayDate') && ~isempty(cfg.dayDate)
                try
                    app.DayofYearDatePicker.Value = datetime(cfg.dayDate, 'InputFormat', 'yyyy-MM-dd');
                catch
                    % Leave the date picker untouched if the stored date can't be parsed
                end
            end
            if isfield(cfg,'timeOfDay'),      app.TimeofdaysecondsEditField.Value = cfg.timeOfDay;        end
            if isfield(cfg,'f107average'),    app.F107averageEditField.Value = cfg.f107average;           end
            if isfield(cfg,'f107daily'),      app.F107dailyEditField.Value = cfg.f107daily;               end
            if isfield(cfg,'apIndex'),        app.APMagneticIndexEditField_2.Value = cfg.apIndex;         end
            if isfield(cfg,'anomalousOxygen'),app.AnomalousOxygenCheckBox.Value = logical(cfg.anomalousOxygen); end
            
            if isfield(cfg,'aoa'), app.AngleofAttackEditField.Value = cfg.aoa; end
            if isfield(cfg,'aos'), app.AngleofSideslipEditField.Value = cfg.aos; end
            if isfield(cfg,'gsiModel') && any(strcmp(cfg.gsiModel, app.GSIModelDropDown.Items))
                app.GSIModelDropDown.Value = cfg.gsiModel;
            end
            if isfield(cfg,'accommodation'), app.EditField.Value = cfg.accommodation; end
            if isfield(cfg,'wallTemp'),      app.WallTermperatureEditField.Value = cfg.wallTemp; end
            if isfield(cfg,'shadow'),         app.ShadowCheckBox.Value = logical(cfg.shadow); end
            if isfield(cfg,'hyperthermal'),   app.HyperthermalCheckBox.Value = logical(cfg.hyperthermal); end
            if isfield(cfg,'rays'),           app.RaysEditField.Value = cfg.rays; end
            if isfield(cfg,'solarRadiation'), app.SolarRadiationPressureCheckBox.Value = logical(cfg.solarRadiation); end
            if isfield(cfg,'specularReflectivity'), app.SpecularReflectivityEditField.Value = cfg.specularReflectivity; end
            if isfield(cfg,'diffuseReflectivity'),  app.DiffuseReflecitvityEditField.Value = cfg.diffuseReflectivity;   end
            if isfield(cfg,'verbose'),    app.VerboseCheckBox.Value = logical(cfg.verbose); end
            if isfield(cfg,'deleteTemp'), app.delCheckBox.Value = logical(cfg.deleteTemp);   end
            if isfield(cfg,'plotField') && any(strcmp(cfg.plotField, app.PlotDropDown.Items))
                app.PlotDropDown.Value = cfg.plotField;
            end
            if isfield(cfg,'lockColorScale')
                app.LockColorScaleCheckBox.Value = logical(cfg.lockColorScale);
                app.LockedCLim = []; % re-arm: capture fresh on next render rather than reuse a stale range
            end
            
            app.logMessage(sprintf('Preset loaded: %s', file));
        end
        
        % Export the plot in the currently active tab as PNG or FIG
        function ExportPlotButtonPushed(app, event)
            activeTab = app.RightTabGroup.SelectedTab;
            
            if activeTab == app.PlotTab2
                % Merged Mesh tab overlays two axes - export the plot thats
                % currently being shown rather than the one that's found
                % first
                if strcmp(app.MeshViewSwitch.Value, 'Final')
                    axSrc = app.UIAxesMeshFinal;
                else
                    axSrc = app.UIAxes2;
                end
            else
                axCandidates = findobj(activeTab, 'Type', 'axes');
                if isempty(axCandidates)
                    uialert(app.UIFigure, 'The active tab has no plot to export.', 'Nothing to Export');
                    return;
                end
                axSrc = axCandidates(1);
            end
            
            defaultName = [matlab.lang.makeValidName(activeTab.Title), '.png'];
            [file, path, filterIdx] = uiputfile({'*.png', 'PNG Image (*.png)'; '*.fig', 'MATLAB Figure (*.fig)'}, ...
                'Export Current Plot As', defaultName);
            if isequal(file, 0)
                return;
            end
            fullPath = fullfile(path, file);
            
            try
                if filterIdx == 2
                    tmpFig = figure('Visible', 'off', 'Color', 'w');
                    tmpAx = copyobj(axSrc, tmpFig);
                    tmpAx.Units = 'normalized';
                    tmpAx.OuterPosition = [0 0 1 1];
                    savefig(tmpFig, fullPath);
                    close(tmpFig);
                else
                    exportgraphics(axSrc, fullPath, 'Resolution', 200);
                end
                app.logMessage(sprintf('Plot exported: %s', fullPath));
            catch ME
                app.logMessage(sprintf('PLOT EXPORT FAILED: %s', ME.message));
            end
        end
        
        % Export the most recent run's results as the raw .mat file or a
        % flattened summary .csv
        function ExportResultsButtonPushed(app, event)
            if isempty(app.LastFileOut) || ~exist(app.LastFileOut, 'file')
                uialert(app.UIFigure, 'No results available yet - run a simulation first.', 'Nothing to Export');
                return;
            end
            
            [file, path, filterIdx] = uiputfile({'*.mat', 'MATLAB Data File (*.mat)'; '*.csv', 'Summary CSV (*.csv)'}, ...
                'Export Results As', 'adbsat_results.mat');
            if isequal(file, 0)
                return;
            end
            fullPath = fullfile(path, file);
            
            try
                if filterIdx == 2
                    if isempty(app.LastSummaryTable) || height(app.LastSummaryTable) == 0
                        error('No results summary is available to export yet.');
                    end
                    writetable(app.LastSummaryTable, fullPath);
                else
                    copyfile(app.LastFileOut, fullPath);
                end
                app.logMessage(sprintf('Results exported: %s', fullPath));
            catch ME
                app.logMessage(sprintf('RESULTS EXPORT FAILED: %s', ME.message));
            end
        end
        
        % Exports every populated plot tab (surface plot, mesh, quality,
        % normals, material ID) to PNG files in one folder, so the user
        % doesn't have to export each tab individually via Export Plot...
        function ExportAllPlotsButtonPushed(app, event)
            folder = uigetdir(pwd, 'Select Folder to Export All Plots Into');
            if isequal(folder, 0)
                return;
            end
            
            if strcmp(app.MeshViewSwitch.Value, 'Final')
                meshAx = app.UIAxesMeshFinal;
            else
                meshAx = app.UIAxes2;
            end
            
            exportList = { ...
                'Surface Distribution', app.UIAxes1; ...
                'Mesh',                 meshAx; ...
                'Mesh Quality',         app.UIAxesQuality; ...
                'Normals',              app.UIAxesNormals; ...
                'Material ID',          app.UIAxesMaterial};
            
            exportedCount = 0;
            for k = 1:size(exportList, 1)
                label = exportList{k, 1};
                ax = exportList{k, 2};
                if isempty(ax.Children)
                    continue; % skip tabs with nothing plotted yet
                end
                fname = [matlab.lang.makeValidName(label), '.png'];
                fullPath = fullfile(folder, fname);
                try
                    exportgraphics(ax, fullPath, 'Resolution', 200);
                    exportedCount = exportedCount + 1;
                catch ME
                    app.logMessage(sprintf('EXPORT FAILED (%s): %s', label, ME.message));
                end
            end
            
            if exportedCount == 0
                app.logMessage('Export All Plots: nothing to export yet.');
            else
                app.logMessage(sprintf('Export All Plots: %d plot(s) exported to %s', exportedCount, folder));
            end
        end
        
        % Rebuilds the Results Summary tab from the output .mat (this is 
        % deliberately generic so it keeps working even if calc_coeff's 
        % output field names change or vary by version, rather than 
        % hard-coding specific coefficient names).
        function updateResultsSummary(app, modOut, fileOut, alt, inc, aoa_deg, aos_deg, inparam, shadow, solar)
            rows = cell(0, 2);
            rows(end+1, :) = {'Model File', app.SelectedPathEditField.Value};
            
            try
                meshS = load(modOut, 'meshdata');
                md = meshS.meshdata;
                rows(end+1, :) = {'Element Count', numel(md.Areas)};
                rows(end+1, :) = {'Total Surface Area (m^2)', sum(md.Areas)};
                rows(end+1, :) = {'Reference Length (m)', md.Lref};
            catch
                rows(end+1, :) = {'Mesh Stats', 'unavailable'};
            end
            
            if shadow, shadowTxt = 'On'; else, shadowTxt = 'Off'; end
            if solar,  solarTxt  = 'On'; else, solarTxt  = 'Off'; end
            if isfield(inparam, 'hyperthermal') && inparam.hyperthermal
                hyperthermalTxt = 'Yes';
            else
                hyperthermalTxt = 'No';
            end
            
            rows(end+1, :) = {'Altitude (km)', alt};
            rows(end+1, :) = {'Inclination (deg)', inc};
            rows(end+1, :) = {'Angle of Attack (deg)', aoa_deg(1)};
            rows(end+1, :) = {'Angle of Sideslip (deg)', aos_deg(1)};
            rows(end+1, :) = {'GSI Model', inparam.gsi_model};
            rows(end+1, :) = {'Hyperthermal Flow Assumed', hyperthermalTxt};
            rows(end+1, :) = {'Accommodation Coeff.', mat2str(inparam.alpha)};
            rows(end+1, :) = {'Wall Temperature (K)', inparam.Tw};
            rows(end+1, :) = {'Shadowing', shadowTxt};
            rows(end+1, :) = {'Solar Radiation Pressure', solarTxt};
            
            try
                s = load(fileOut);
                fn = fieldnames(s);
                for i = 1:numel(fn)
                    v = s.(fn{i});
                    if isnumeric(v) && isscalar(v)
                        rows(end+1, :) = {['Result: ', fn{i}], v}; %#ok<AGROW>
                    end
                end
            catch

            end
            
            app.ResultsTable.Data = rows;
            valueStrings = cellfun(@(v) string(v), rows(:,2), 'UniformOutput', true);
            app.LastSummaryTable = table(string(rows(:,1)), valueStrings, 'VariableNames', {'Parameter', 'Value'});
        end

        % Triggered when Cancel Button is hit during a run
        function CancelButtonPushed(app, event)
            app.CancelRequested = true;
            app.logMessage('Cancel requested... awaiting safe checkpoint to terminate.');
            app.CancelButton.Enable = 'off';
        end

        % Main solver execution workflow invocation logic.
        %
        % NOTE: This button now runs everything from CHECKPOINT 2 onward
        % only. The geometry import / mesh generation stage (previously 
        % Checkpoint 1 -> 2) must be run first via the "Import & Generate 
        % Mesh" button on the Model Geometry tab
        function RunButtonPushed(app, event)
            
            % Enforce button states and flags
            app.RunButton.Enable = 'off';
            app.CancelButton.Enable = 'on';
            app.CancelRequested = false;
            app.ProgressBar.Value = 0;

            app.clearLog();
            app.logMessage('>> Initializing pipeline execution loops...');
            

            app.ProgressBar.Value = 5; drawnow;

            % A processed mesh must already exist - generated via the
            % "Import & Generate Mesh" button on the Model Geometry tab -
            % since this button no longer performs that step itself
            if isempty(app.LastModOut) || ~exist(app.LastModOut, 'file')
                app.logMessage('ABORT: No processed mesh available. Click "Import & Generate Mesh" on the Model Geometry tab first.');
                uialert(app.UIFigure, 'Please import & generate a mesh (Model Geometry tab) before running the simulation.', 'No Mesh Available');
                app.RunButton.Enable = 'on';
                app.CancelButton.Enable = 'off';
                return;
            end
            
            modOut = app.LastModOut;
            [~, modName] = fileparts(modOut);
            
            % Only clear the aerodynamic result plot - the mesh/quality/
            % normals/material plots from the last import are still valid
            % and are left alone
            app.clearAerodynamicPlot();
            
            try
                ADBSat_path = ADBSat_dynpath;
                resOut = fullfile(ADBSat_path, 'inou', 'results', modName);
            catch
                resOut = fullfile(fileparts(modOut), [modName, '_results']);
            end
            
            app.logMessage(sprintf('Processing: %s', upper(modName)));
            app.ProgressBar.Value = 10; drawnow;
            
            % Early cancel check before parameter validation
            if app.CancelRequested
                app.logMessage('Simulation safely aborted by user.');
                app.RunButton.Enable = 'on';
                app.CancelButton.Enable = 'off';
                return;
            end
            
            % Collect UI Control parameter boundaries
            alt = app.AltitudekmEditField.Value;
            inc = app.InclinationdegEditField.Value;
            
            % -----------------------------------------------------------
            % HARD INPUT VALIDATION - abort the run rather than let these
            % fail deep inside the solver pipeline with a confusing
            % low-level error
            % -----------------------------------------------------------
            validationErrors = {};
            if isnan(alt) || alt <= 0
                validationErrors{end+1} = 'Altitude (km) must be a positive number.';
            end
            
            if ~isempty(validationErrors)
                for k = 1:numel(validationErrors)
                    app.logMessage(sprintf('ABORT: %s', validationErrors{k}));
                end
                uialert(app.UIFigure, strjoin(validationErrors, newline), 'Invalid Input');
                app.RunButton.Enable = 'on';
                app.CancelButton.Enable = 'off';
                return;
            end
            
            dpVal = app.DayofYearDatePicker.Value;
            if ~isempty(dpVal)
                dayOfYear = day(datetime(dpVal), 'dayofyear');
            else
                dayOfYear = 1; 
            end
            
            ap_val = str2num(app.APMagneticIndexEditField_2.Value);
            if isempty(ap_val) || numel(ap_val) < 7
                app.logMessage('WARNING: AP Magnetic Index needs 7 values - defaulting to all zeros.');
                ap_val = zeros(1, 7);
            end
            
            env = [alt*1e3, inc/2, app.LongitudedegEditField.Value, ...
                   dayOfYear, app.TimeofdaysecondsEditField.Value, ...
                   app.F107dailyEditField.Value, app.F107averageEditField.Value, ...
                   ap_val, app.AnomalousOxygenCheckBox.Value];
               
            aoa_deg = app.AngleofAttackEditField.Value;
            aos_deg = app.AngleofSideslipEditField.Value;
            shadow = app.ShadowCheckBox.Value;
            inparam.gsi_model = app.GSIModelDropDown.Value;
            inparam.hyperthermal = app.HyperthermalCheckBox.Value;
            inparam.rays = app.RaysEditField.Value;
            
            alpha_val = str2num(app.EditField.Value);
            if isempty(alpha_val)
                app.logMessage('WARNING: Accommodation Coefficient could not be parsed - defaulting to 1.');
                alpha_val = 1;
            elseif any(alpha_val < 0 | alpha_val > 1)
                app.logMessage('WARNING: Accommodation Coefficient(s) outside the expected [0,1] range - proceeding anyway.');
            end
            inparam.alpha = alpha_val;
            
            tw_val = str2double(app.WallTermperatureEditField.Value);
            if isnan(tw_val) || tw_val <= 0
                app.logMessage('WARNING: Wall Temperature invalid - defaulting to 300 K.');
                tw_val = 300;
            end
            inparam.Tw = tw_val;
            
            solar = app.SolarRadiationPressureCheckBox.Value;
            inparam.sol_cR = app.SpecularReflectivityEditField.Value;
            inparam.sol_cD = app.DiffuseReflecitvityEditField.Value;
            
            verb = app.VerboseCheckBox.Value;
            del = app.delCheckBox.Value;
            
            app.ProgressBar.Value = 50; drawnow;
            
            % CHECKPOINT 2
            if app.CancelRequested
                app.logMessage('Simulation safely aborted by user before aerodynamic computations.');
                app.RunButton.Enable = 'on';
                app.CancelButton.Enable = 'off';
                return;
            end
            
            % Calculate aerodynamics coefficients
            app.logMessage('Computing aerodynamic flux and interactions...');
            try
                inparam = environment(inparam, env(1), env(2), env(3), env(4), env(5), env(6), env(7), env(8:14), env(15));
                fileOut = calc_coeff(modOut, resOut, deg2rad(aoa_deg), deg2rad(aos_deg), inparam, shadow, solar, 1, 0); 
                app.logMessage('Aerodynamic evaluation routines mapped.');
            catch ME
                app.logMessage(sprintf('SOLVER FAULT: %s', ME.message));
                app.RunButton.Enable = 'on';
                app.CancelButton.Enable = 'off';
                return;
            end
            
            app.LastFileOut = fileOut;
            try
                app.updateResultsSummary(modOut, fileOut, alt, inc, aoa_deg, aos_deg, inparam, shadow, solar);
            catch ME
                app.logMessage(sprintf('Could not build results summary: %s', ME.message));
            end
            
            app.ProgressBar.Value = 85; drawnow;

            % CHECKPOINT 3
            if app.CancelRequested
                app.logMessage('Simulation safely aborted by user before plotting outputs.');
                app.RunButton.Enable = 'on';
                app.CancelButton.Enable = 'off';
                return;
            end
            
            % -----------------------------------------------------------
            % 3. PLOT SURFACE RESPONSE RESULTS
            % -----------------------------------------------------------
            if verb && ~del
                app.logMessage('Mapping data arrays onto 3D graphics axes layers...');
                
                app.renderSurfacePlot(fileOut, modOut, aoa_deg(1), aos_deg(1), app.PlotDropDown.Value);
                
                app.RightTabGroup.SelectedTab = app.PlotTab1;
            end
            
            app.ProgressBar.Value = 100; drawnow;
            app.logMessage('>> Pipeline evaluation complete.');
            
            % Reset UI Elements once complete
            app.RunButton.Enable = 'on';
            app.CancelButton.Enable = 'off';
        end
        
        % -----------------------------------------------------------
        % Shared geometry import + mesh generation stage (everything up
        % to CHECKPOINT 2)
        % -----------------------------------------------------------
        function [modOut, success] = runGeometryImport(app, modIn, ADBSat_path, verb, meshParam, axs)
            modOut = '';
            success = false;
            
            targetDir = fullfile(ADBSat_path, 'inou', 'models');
            if ~exist(targetDir, 'dir')
                app.logMessage('Directory "inou/models" not found. Creating it now...');
                mkdir(targetDir);
            end
            
            app.logMessage('Parsing mesh via ADBSat geometry engine...');
            try
                modOut = ADBSatImport_app(modIn, targetDir, verb, meshParam, axs);
                app.logMessage('Mesh conversion and geometry tracking processing finalized.');
                app.applyMeshViewToggle();
                success = true;
            catch ME
                app.logMessage(sprintf('GEOMETRY PARSING FAILURE: %s', ME.message));
                success = false;
            end
        end
        
        % Standalone entry point on the Model Geometry tab: runs the
        % import + mesh-generation stage only (up to CHECKPOINT 2),
        % without the aerodynamic solver or result plotting stages.
        % Useful for previewing/validating a mesh before committing to a
        % full simulation run
        function ProcessGeometryButtonPushed(app, event)
            app.ProcessGeometryButton.Enable = 'off';
            app.RunButton.Enable = 'off';
            
            app.clearLog();
            app.logMessage('>> Running geometry import & mesh generation only...');
            app.ProgressBar.Value = 0; drawnow;
            
            if isempty(app.UserSelectedFullFile) || ~exist(app.UserSelectedFullFile, 'file')
                app.logMessage('ABORT: No valid input geometry model file linked! Click Browse first.');
                uialert(app.UIFigure, 'Please select an input geometry model file before attempting to process it.', 'Missing Model File');
                app.ProcessGeometryButton.Enable = 'on';
                app.RunButton.Enable = 'on';
                return;
            end
            
            app.clearAllPlots();
            modIn = app.UserSelectedFullFile;
            
            elementSizeVal = str2double(app.ElementSizemEditField.Value);
            subdivisionsVal = app.NumberofSubdivisionsEditField.Value;
            
            validationErrors = {};
            if isnan(elementSizeVal) || elementSizeVal <= 0
                validationErrors{end+1} = 'Element Size (m) must be a positive number.'; %#ok<AGROW>
            end
            if isnan(subdivisionsVal) || subdivisionsVal < 0 || mod(subdivisionsVal, 1) ~= 0
                validationErrors{end+1} = 'Number of Subdivisions must be a non-negative whole number.'; %#ok<AGROW>
            end
            
            if ~isempty(validationErrors)
                for k = 1:numel(validationErrors)
                    app.logMessage(sprintf('ABORT: %s', validationErrors{k}));
                end
                uialert(app.UIFigure, strjoin(validationErrors, newline), 'Invalid Input');
                app.ProcessGeometryButton.Enable = 'on';
                app.RunButton.Enable = 'on';
                return;
            end
            
            meshParam.subdivisions = subdivisionsVal;
            meshParam.elementSize = elementSizeVal;
            meshParam.translation = [app.inXEditField.Value, app.inYEditField.Value, app.inZEditField.Value];
            meshParam.rotationAngles = [app.PitchEditField.Value, app.RollEditField.Value, app.YawEditField.Value];
            meshParam.centering = app.AutoCenteringCheckBox.Value;
            
            verb = app.VerboseCheckBox.Value;
            
            app.ProgressBar.Value = 20; drawnow;
            
            try
                ADBSat_path = ADBSat_dynpath;
            catch
                ADBSat_path = pwd;
                app.logMessage('WARNING: Could not resolve ADBSat root via ADBSat_dynpath - falling back to current working directory.');
            end
            
            axs.mesh1    = app.UIAxes2;
            axs.mesh2    = app.UIAxesMeshFinal;
            axs.quality  = app.UIAxesQuality;
            axs.normals  = app.UIAxesNormals;
            axs.material = app.UIAxesMaterial;
            
            [modOut, importOk] = app.runGeometryImport(modIn, ADBSat_path, verb, meshParam, axs);
            
            app.ProgressBar.Value = 100; drawnow;
            
            if importOk
                app.LastModOut = modOut;
                app.logMessage('>> Geometry import & mesh generation complete (Checkpoint 2 reached).');
            else
                app.logMessage('>> Geometry import & mesh generation failed.');
            end
            
            app.ProcessGeometryButton.Enable = 'on';
            app.RunButton.Enable = 'on';
        end
        
        % Dynamic Layout
        function updateAppLayout(app, event)
            currWidth = app.UIFigure.Position(3);
            if(currWidth <= app.onePanelWidth)
                app.MainGridLayout.RowHeight = {660, 400, 150};
                app.MainGridLayout.ColumnWidth = {'1x'};
                
                app.LeftPanel.Layout.Row = 1;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
                app.BottomPanel.Layout.Row = 3;
                app.BottomPanel.Layout.Column = 1;
            else
                app.MainGridLayout.RowHeight = {'1x', 140};
                app.MainGridLayout.ColumnWidth = {345, '1x'};
                
                app.LeftPanel.Layout.Row = [1 2];
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
                app.BottomPanel.Layout.Row = 2;
                app.BottomPanel.Layout.Column = 2;
            end
        end
    end
    
    methods (Access = private)
        function createComponents(app)
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 1050 780];
            app.UIFigure.Name = 'ADBSat Spacecraft Aerodynamics Toolkit App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);
            
            app.MainGridLayout = uigridlayout(app.UIFigure);
            app.MainGridLayout.ColumnWidth = {345, '1x'};
            app.MainGridLayout.RowHeight = {'1x', 140};
            app.MainGridLayout.ColumnSpacing = 6;
            app.MainGridLayout.RowSpacing = 6;
            app.MainGridLayout.Padding = [6 6 6 6];
            app.MainGridLayout.Scrollable = 'on';
            
            % Left control options panel config
            app.LeftPanel = uipanel(app.MainGridLayout);
            app.LeftPanel.Layout.Row = [1 2];
            app.LeftPanel.Layout.Column = 1;
            
            % Flexibly fill left panel housing boundaries
            app.SavePresetButton = uibutton(app.LeftPanel, 'push');
            app.SavePresetButton.Position = [4 738 163 26];
            app.SavePresetButton.Text = 'Save Preset...';
            app.SavePresetButton.ButtonPushedFcn = createCallbackFcn(app, @SavePresetButtonPushed, true);
            
            app.LoadPresetButton = uibutton(app.LeftPanel, 'push');
            app.LoadPresetButton.Position = [172 738 163 26];
            app.LoadPresetButton.Text = 'Load Preset...';
            app.LoadPresetButton.ButtonPushedFcn = createCallbackFcn(app, @LoadPresetButtonPushed, true);
            
            app.TabGroup = uitabgroup(app.LeftPanel);
            app.TabGroup.Position = [4 4 335 730];
            
            % -----------------------------------------------------------
            % MODEL SELECTION TAB (COMPACT TOP ALIGNED LAYOUT)
            % -----------------------------------------------------------
            app.ModelTab = uitab(app.TabGroup);
            app.ModelTab.Title = 'Model Geometry';
            
            app.SelectedPathLabel = uilabel(app.ModelTab);
            app.SelectedPathLabel.Position = [10 670 200 22]; % Shifted down
            app.SelectedPathLabel.Text = 'Selected Geometry File:';
            app.SelectedPathLabel.FontWeight = 'bold';
            
            app.SelectedPathEditField = uieditfield(app.ModelTab, 'text');
            app.SelectedPathEditField.Position = [10 643 230 24];
            app.SelectedPathEditField.Editable = 'off';
            app.SelectedPathEditField.Placeholder = 'No model file selected...';
            
            app.BrowseButton = uibutton(app.ModelTab, 'push');
            app.BrowseButton.Position = [246 643 75 24];
            app.BrowseButton.Text = 'Browse...';
            app.BrowseButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseButtonPushed, true);
            
            app.EditField_3Label = uilabel(app.ModelTab);
            app.EditField_3Label.HorizontalAlignment = 'right';
            app.EditField_3Label.Position = [10 600 140 22];
            app.EditField_3Label.Text = 'Number of Subdivisions:';
            
            app.NumberofSubdivisionsEditField = uieditfield(app.ModelTab, 'numeric');
            app.NumberofSubdivisionsEditField.Position = [160 600 160 22];
            app.NumberofSubdivisionsEditField.Value = 0;
            
            app.EditField_3Label_2 = uilabel(app.ModelTab);
            app.EditField_3Label_2.HorizontalAlignment = 'right';
            app.EditField_3Label_2.Position = [10 560 140 22];
            app.EditField_3Label_2.Text = 'Element Size (m):';
            
            app.ElementSizemEditField = uieditfield(app.ModelTab, 'text');
            app.ElementSizemEditField.Position = [160 560 160 22];
            app.ElementSizemEditField.Value = '0.002';
            
            app.AutoCenteringCheckBox = uicheckbox(app.ModelTab);
            app.AutoCenteringCheckBox.Text = 'Auto Centering Mesh Object';
            app.AutoCenteringCheckBox.Position = [10 525 250 22];
            app.AutoCenteringCheckBox.Value = true;
            
            app.TranslationLabel = uilabel(app.ModelTab);
            app.TranslationLabel.Position = [20 480 100 22];
            app.TranslationLabel.Text = 'Translation';
            app.TranslationLabel.FontWeight = 'bold';
            
            app.EditField_3Label_3 = uilabel(app.ModelTab);
            app.EditField_3Label_3.HorizontalAlignment = 'right';
            app.EditField_3Label_3.Position = [10 450 35 22];
            app.EditField_3Label_3.Text = 'in X:';
            app.inXEditField = uieditfield(app.ModelTab, 'numeric');
            app.inXEditField.Position = [50 450 50 22];
            
            app.inYLabel = uilabel(app.ModelTab);
            app.inYLabel.HorizontalAlignment = 'right';
            app.inYLabel.Position = [115 450 35 22];
            app.inYLabel.Text = 'in Y:';
            app.inYEditField = uieditfield(app.ModelTab, 'numeric');
            app.inYEditField.Position = [155 450 50 22];
            
            app.inZLabel = uilabel(app.ModelTab);
            app.inZLabel.HorizontalAlignment = 'right';
            app.inZLabel.Position = [220 450 35 22];
            app.inZLabel.Text = 'in Z:';
            app.inZEditField = uieditfield(app.ModelTab, 'numeric');
            app.inZEditField.Position = [260 450 60 22];
            
            app.RotationLabel = uilabel(app.ModelTab);
            app.RotationLabel.Position = [20 405 100 22];
            app.RotationLabel.Text = 'Rotation';
            app.RotationLabel.FontWeight = 'bold';
            
            app.PitchEditFieldLabel = uilabel(app.ModelTab);
            app.PitchEditFieldLabel.HorizontalAlignment = 'right';
            app.PitchEditFieldLabel.Position = [10 375 38 22];
            app.PitchEditFieldLabel.Text = 'Pitch:';
            app.PitchEditField = uieditfield(app.ModelTab, 'numeric');
            app.PitchEditField.Position = [52 375 48 22];
            
            app.RollLabel = uilabel(app.ModelTab);
            app.RollLabel.HorizontalAlignment = 'right';
            app.RollLabel.Position = [115 375 35 22];
            app.RollLabel.Text = 'Roll:';
            app.RollEditField = uieditfield(app.ModelTab, 'numeric');
            app.RollEditField.Position = [155 375 50 22];
            
            app.YawEditFieldLabel = uilabel(app.ModelTab);
            app.YawEditFieldLabel.HorizontalAlignment = 'right';
            app.YawEditFieldLabel.Position = [220 375 35 22];
            app.YawEditFieldLabel.Text = 'Yaw:';
            app.YawEditField = uieditfield(app.ModelTab, 'numeric');
            app.YawEditField.Position = [260 375 60 22];
            
            % Runs import + mesh generation only (through Checkpoint 2),
            % without running the full aerodynamic solver / plotting pipeline
            app.ProcessGeometryButton = uibutton(app.ModelTab, 'push');
            app.ProcessGeometryButton.Position = [10 320 310 32];
            app.ProcessGeometryButton.Text = 'Import & Generate Mesh';
            app.ProcessGeometryButton.FontWeight = 'bold';
            app.ProcessGeometryButton.BackgroundColor = [0.2 0.4 0.7];
            app.ProcessGeometryButton.FontColor = [1 1 1];
            app.ProcessGeometryButton.ButtonPushedFcn = createCallbackFcn(app, @ProcessGeometryButtonPushed, true);
            
            % -----------------------------------------------------------
            % ENVIRONMENT TAB 
            % -----------------------------------------------------------
            app.EnvironmentTab = uitab(app.TabGroup);
            app.EnvironmentTab.Title = 'Environment';
            
            app.AltitudekmEditFieldLabel = uilabel(app.EnvironmentTab);
            app.AltitudekmEditFieldLabel.HorizontalAlignment = 'right';
            app.AltitudekmEditFieldLabel.Position = [10 670 120 22];
            app.AltitudekmEditFieldLabel.Text = 'Altitude (km):';
            app.AltitudekmEditField = uieditfield(app.EnvironmentTab, 'numeric');
            app.AltitudekmEditField.Position = [140 670 180 22];
            app.AltitudekmEditField.Value = 400;
            
            app.InclinationdegEditFieldLabel = uilabel(app.EnvironmentTab);
            app.InclinationdegEditFieldLabel.HorizontalAlignment = 'right';
            app.InclinationdegEditFieldLabel.Position = [10 630 120 22];
            app.InclinationdegEditFieldLabel.Text = 'Inclination (deg):';
            app.InclinationdegEditField = uieditfield(app.EnvironmentTab, 'numeric');
            app.InclinationdegEditField.Position = [140 630 180 22];
            
            app.LongitudedegEditFieldLabel = uilabel(app.EnvironmentTab);
            app.LongitudedegEditFieldLabel.HorizontalAlignment = 'right';
            app.LongitudedegEditFieldLabel.Position = [10 590 120 22];
            app.LongitudedegEditFieldLabel.Text = 'Longitude (deg):';
            app.LongitudedegEditField = uieditfield(app.EnvironmentTab, 'numeric');
            app.LongitudedegEditField.Position = [140 590 180 22];
            
            app.DayofYearDatePickerLabel = uilabel(app.EnvironmentTab);
            app.DayofYearDatePickerLabel.HorizontalAlignment = 'right';
            app.DayofYearDatePickerLabel.Position = [10 550 120 22];
            app.DayofYearDatePickerLabel.Text = 'Day of Year:';
            app.DayofYearDatePicker = uidatepicker(app.EnvironmentTab);
            app.DayofYearDatePicker.Position = [140 550 180 22];
            
            app.TimeofdaysecondsEditFieldLabel = uilabel(app.EnvironmentTab);
            app.TimeofdaysecondsEditFieldLabel.HorizontalAlignment = 'right';
            app.TimeofdaysecondsEditFieldLabel.Position = [5 510 125 22];
            app.TimeofdaysecondsEditFieldLabel.Text = 'Time of Day (sec):';
            app.TimeofdaysecondsEditField = uieditfield(app.EnvironmentTab, 'numeric');
            app.TimeofdaysecondsEditField.Position = [140 510 180 22];
            
            app.F107averageEditFieldLabel = uilabel(app.EnvironmentTab);
            app.F107averageEditFieldLabel.HorizontalAlignment = 'right';
            app.F107averageEditFieldLabel.Position = [10 470 120 22];
            app.F107averageEditFieldLabel.Text = 'F10.7 Average:';
            app.F107averageEditField = uieditfield(app.EnvironmentTab, 'numeric');
            app.F107averageEditField.Position = [140 470 180 22];
            
            app.F107dailyEditFieldLabel = uilabel(app.EnvironmentTab);
            app.F107dailyEditFieldLabel.HorizontalAlignment = 'right';
            app.F107dailyEditFieldLabel.Position = [10 430 120 22];
            app.F107dailyEditFieldLabel.Text = 'F10.7 Daily:';
            app.F107dailyEditField = uieditfield(app.EnvironmentTab, 'numeric');
            app.F107dailyEditField.Position = [140 430 180 22];
            
            app.EditFieldLabel = uilabel(app.EnvironmentTab);
            app.EditFieldLabel.HorizontalAlignment = 'right';
            app.EditFieldLabel.Position = [10 390 120 22];
            app.EditFieldLabel.Text = 'AP Magnetic Index:';
            app.APMagneticIndexEditField_2 = uieditfield(app.EnvironmentTab, 'text');
            app.APMagneticIndexEditField_2.Position = [140 390 180 22];
            app.APMagneticIndexEditField_2.Value = '0 0 0 0 0 0 0';
            
            app.AnomalousOxygenCheckBox = uicheckbox(app.EnvironmentTab);
            app.AnomalousOxygenCheckBox.Text = 'Enable Anomalous Oxygen Model';
            app.AnomalousOxygenCheckBox.Position = [10 355 280 22];
            
            % -----------------------------------------------------------
            % SIMULATION CONFIG TAB
            % -----------------------------------------------------------
            app.SimulationTab = uitab(app.TabGroup);
            app.SimulationTab.Title = 'Simulation';
            
            % Container to hold controls, leaving space at the bottom
            % This keeps controls from overlapping the progress bar
            simControlPanel = uipanel(app.SimulationTab);
            simControlPanel.Position = [0 40 335 680]; 
            simControlPanel.BorderType = 'none';

            % Progress Bar - Placed explicitly in the bottom 30px of the Tab
            app.ProgressBar = uigauge(app.SimulationTab, 'linear');
            app.ProgressBar.Position = [10 10 315 25]; % [Left Bottom Width Height]
            app.ProgressBar.Limits = [0 100];
            app.ProgressBar.Value = 0;
            app.ProgressBar.MajorTicks = [];
            
            app.AngleofAttackEditFieldLabel = uilabel(simControlPanel);
            app.AngleofAttackEditFieldLabel.HorizontalAlignment = 'right';
            app.AngleofAttackEditFieldLabel.Position = [10 630 110 22];
            app.AngleofAttackEditFieldLabel.Text = 'Angle of Attack:';
            app.AngleofAttackEditField = uieditfield(simControlPanel, 'numeric');
            app.AngleofAttackEditField.Limits = [0 360];
            app.AngleofAttackEditField.Position = [130 630 190 22];
            app.AngleofAttackEditField.Value = 0;
            
            app.AngleofSideslipEditFieldLabel = uilabel(simControlPanel);
            app.AngleofSideslipEditFieldLabel.HorizontalAlignment = 'right';
            app.AngleofSideslipEditFieldLabel.Position = [10 590 110 22];
            app.AngleofSideslipEditFieldLabel.Text = 'Angle of Sideslip:';
            app.AngleofSideslipEditField = uieditfield(simControlPanel, 'numeric');
            app.AngleofSideslipEditField.Limits = [0 360];
            app.AngleofSideslipEditField.Position = [130 590 190 22];
            app.AngleofSideslipEditField.Value = 0;
            
            app.GSIModelDropDownLabel = uilabel(simControlPanel);
            app.GSIModelDropDownLabel.HorizontalAlignment = 'right';
            app.GSIModelDropDownLabel.Position = [10 550 110 22];
            app.GSIModelDropDownLabel.Text = 'GSI Force Model:';
            app.GSIModelDropDown = uidropdown(simControlPanel);
            app.GSIModelDropDown.Items = {'DRIA', 'cook', 'schaaf', 'sentman', 'CLL', 'maxwell', 'newton', 'storchHyp'};
            app.GSIModelDropDown.Position = [130 550 190 22];
            app.GSIModelDropDown.Value = 'DRIA';
            
            app.AccomodationCoefficientsLabel = uilabel(simControlPanel);
            app.AccomodationCoefficientsLabel.HorizontalAlignment = 'right';
            app.AccomodationCoefficientsLabel.Position = [10 510 110 22];
            app.AccomodationCoefficientsLabel.Text = 'Accomodation Coeff:';
            app.EditField = uieditfield(simControlPanel, 'text');
            app.EditField.Position = [130 510 190 22];
            app.EditField.Value = '1';
            
            app.WallTermperatureKelvinLabel = uilabel(simControlPanel);
            app.WallTermperatureKelvinLabel.HorizontalAlignment = 'right';
            app.WallTermperatureKelvinLabel.Position = [10 470 110 22];
            app.WallTermperatureKelvinLabel.Text = 'Wall Temp (K):';
            app.WallTermperatureEditField = uieditfield(simControlPanel, 'text');
            app.WallTermperatureEditField.Position = [130 470 190 22];
            app.WallTermperatureEditField.Value = '300';
            
            app.ShadowCheckBox = uicheckbox(simControlPanel);
            app.ShadowCheckBox.Text = 'Shadowing';
            app.ShadowCheckBox.Position = [130 435 90 22];
            
            app.HyperthermalCheckBox = uicheckbox(simControlPanel);
            app.HyperthermalCheckBox.Text = 'Hyperthermal';
            app.HyperthermalCheckBox.Position = [230 435 90 22];
            
            app.RaysEditFieldLabel = uilabel(simControlPanel);
            app.RaysEditFieldLabel.HorizontalAlignment = 'right';
            app.RaysEditFieldLabel.Position = [10 400 110 22];
            app.RaysEditFieldLabel.Text = 'Number of Rays:';
            app.RaysEditField = uieditfield(simControlPanel, 'numeric');
            app.RaysEditField.Limits = [1 Inf];
            app.RaysEditField.RoundFractionalValues = 'on';
            app.RaysEditField.Position = [130 400 190 22];
            app.RaysEditField.Value = 40;
            
            app.SolarRadiationPressureCheckBox = uicheckbox(simControlPanel);
            app.SolarRadiationPressureCheckBox.Text = 'Enable Solar Radiation Pressure';
            app.SolarRadiationPressureCheckBox.Position = [130 365 190 22];
            
            app.EditField_2Label = uilabel(simControlPanel);
            app.EditField_2Label.HorizontalAlignment = 'right';
            app.EditField_2Label.Position = [10 330 110 22];
            app.EditField_2Label.Text = 'Specular Reflectivity:';
            app.SpecularReflectivityEditField = uieditfield(simControlPanel, 'numeric');
            app.SpecularReflectivityEditField.Position = [130 330 190 22];
            
            app.EditField_2Label_2 = uilabel(simControlPanel);
            app.EditField_2Label_2.HorizontalAlignment = 'right';
            app.EditField_2Label_2.Position = [10 290 110 22];
            app.EditField_2Label_2.Text = 'Diffuse Reflectivity:';
            app.DiffuseReflecitvityEditField = uieditfield(simControlPanel, 'numeric');
            app.DiffuseReflecitvityEditField.Position = [130 290 190 22];
            
            app.VerboseCheckBox = uicheckbox(simControlPanel);
            app.VerboseCheckBox.Text = 'Verbose Logs';
            app.VerboseCheckBox.Position = [130 255 100 22];
            app.VerboseCheckBox.Value = true;
            
            app.delCheckBox = uicheckbox(simControlPanel);
            app.delCheckBox.Text = 'Delete Temp';
            app.delCheckBox.Position = [240 255 90 22];
            
            app.PlotDropDownLabel = uilabel(simControlPanel);
            app.PlotDropDownLabel.HorizontalAlignment = 'right';
            app.PlotDropDownLabel.Position = [10 210 110 22];
            app.PlotDropDownLabel.Text = 'Active Plot Field:';
            app.PlotDropDown = uidropdown(simControlPanel);
            app.PlotDropDown.Items = {'cl', 'cd', 'cp', 'ctau', 'delta', 'shadow'};
            app.PlotDropDown.Position = [130 210 190 22];
            app.PlotDropDown.Value = 'cl';
            app.PlotDropDown.ValueChangedFcn = createCallbackFcn(app, @PlotDropDownValueChanged, true);
            
            app.LockColorScaleCheckBox = uicheckbox(simControlPanel);
            app.LockColorScaleCheckBox.Text = 'Lock Color Scale';
            app.LockColorScaleCheckBox.Position = [130 184 190 22];
            app.LockColorScaleCheckBox.ValueChangedFcn = createCallbackFcn(app, @LockColorScaleCheckBoxValueChanged, true);
            
            % RESIZED RUN BUTTON
            app.RunButton = uibutton(simControlPanel, 'push');
            app.RunButton.Position = [130 150 90 32];
            app.RunButton.Text = 'Run Sim';
            app.RunButton.FontWeight = 'bold';
            app.RunButton.BackgroundColor = [0.2 0.6 0.2];
            app.RunButton.FontColor = [1 1 1];
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);

            % NEW CANCEL BUTTON
            app.CancelButton = uibutton(simControlPanel, 'push');
            app.CancelButton.Position = [230 150 90 32];
            app.CancelButton.Text = 'Cancel';
            app.CancelButton.FontWeight = 'bold';
            app.CancelButton.BackgroundColor = [0.8 0.2 0.2];
            app.CancelButton.FontColor = [1 1 1];
            app.CancelButton.Enable = 'off';
            app.CancelButton.ButtonPushedFcn = createCallbackFcn(app, @CancelButtonPushed, true);
            
            app.ClearPlotsButton = uibutton(simControlPanel, 'push');
            app.ClearPlotsButton.Position = [130 110 190 28];
            app.ClearPlotsButton.Text = 'Clear Plots';
            app.ClearPlotsButton.ButtonPushedFcn = createCallbackFcn(app, @ClearPlotsButtonPushed, true);
            
            app.ExportPlotButton = uibutton(simControlPanel, 'push');
            app.ExportPlotButton.Position = [10 70 150 28];
            app.ExportPlotButton.Text = 'Export Plot...';
            app.ExportPlotButton.ButtonPushedFcn = createCallbackFcn(app, @ExportPlotButtonPushed, true);
            
            app.ExportResultsButton = uibutton(simControlPanel, 'push');
            app.ExportResultsButton.Position = [170 70 150 28];
            app.ExportResultsButton.Text = 'Export Results...';
            app.ExportResultsButton.ButtonPushedFcn = createCallbackFcn(app, @ExportResultsButtonPushed, true);
            
            app.ExportAllPlotsButton = uibutton(simControlPanel, 'push');
            app.ExportAllPlotsButton.Position = [10 30 310 28];
            app.ExportAllPlotsButton.Text = 'Export All Plots...';
            app.ExportAllPlotsButton.ButtonPushedFcn = createCallbackFcn(app, @ExportAllPlotsButtonPushed, true);
            
            % -----------------------------------------------------------
            % RIGHT TABS PANEL
            % -----------------------------------------------------------
            app.RightPanel = uipanel(app.MainGridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;
            
            rightGrid = uigridlayout(app.RightPanel, [1, 1], 'Padding', [0 0 0 0]);
            
            app.RightTabGroup = uitabgroup(rightGrid);
            app.RightTabGroup.Layout.Row = 1;
            app.RightTabGroup.Layout.Column = 1;
            
            app.PlotTab1 = uitab(app.RightTabGroup);
            app.PlotTab1.Title = '3D Surface Plot';
            tab1Grid = uigridlayout(app.PlotTab1, [1, 1], 'Padding', [5 5 5 5]);
            app.UIAxes1 = uiaxes(tab1Grid);
            app.UIAxes1.Layout.Row = 1;
            app.UIAxes1.Layout.Column = 1;
            title(app.UIAxes1, 'Surface Distribution');
            
            app.PlotTab2 = uitab(app.RightTabGroup);
            app.PlotTab2.Title = 'Mesh';
            tab2Grid = uigridlayout(app.PlotTab2, [2, 1], 'Padding', [8 8 8 4]);
            tab2Grid.RowHeight = {45, '1x'};
            tab2Grid.RowSpacing = 4;
            
            app.MeshViewSwitch = uiswitch(tab2Grid, 'slider');
            app.MeshViewSwitch.Items = {'Original', 'Final'};
            app.MeshViewSwitch.Value = 'Original';
            app.MeshViewSwitch.Layout.Row = 1;
            app.MeshViewSwitch.Layout.Column = 1;
            app.MeshViewSwitch.ValueChangedFcn = createCallbackFcn(app, @MeshViewSwitchValueChanged, true);
            
            % Both mesh axes occupy the same grid cell and are overlaid;
            % only one is ever visible at a time (toggled via the switch)
            app.UIAxes2 = uiaxes(tab2Grid);
            app.UIAxes2.Layout.Row = 2;
            app.UIAxes2.Layout.Column = 1;
            
            app.UIAxesMeshFinal = uiaxes(tab2Grid);
            app.UIAxesMeshFinal.Layout.Row = 2;
            app.UIAxesMeshFinal.Layout.Column = 1;
            app.UIAxesMeshFinal.Visible = 'off';
            
            app.PlotTab4 = uitab(app.RightTabGroup);
            app.PlotTab4.Title = 'Mesh Quality';
            tab4Grid = uigridlayout(app.PlotTab4, [1, 1], 'Padding', [5 5 5 5]);
            app.UIAxesQuality = uiaxes(tab4Grid);
            app.UIAxesQuality.Layout.Row = 1;
            app.UIAxesQuality.Layout.Column = 1;
            
            app.PlotTab5 = uitab(app.RightTabGroup);
            app.PlotTab5.Title = 'Normals';
            tab5Grid = uigridlayout(app.PlotTab5, [1, 1], 'Padding', [5 5 5 5]);
            app.UIAxesNormals = uiaxes(tab5Grid);
            app.UIAxesNormals.Layout.Row = 1;
            app.UIAxesNormals.Layout.Column = 1;
            
            app.PlotTab6 = uitab(app.RightTabGroup);
            app.PlotTab6.Title = 'Material ID';
            tab6Grid = uigridlayout(app.PlotTab6, [1, 1], 'Padding', [5 5 5 5]);
            app.UIAxesMaterial = uiaxes(tab6Grid);
            app.UIAxesMaterial.Layout.Row = 1;
            app.UIAxesMaterial.Layout.Column = 1;
            
            app.PlotTab7 = uitab(app.RightTabGroup);
            app.PlotTab7.Title = 'Results Summary';
            tab7Grid = uigridlayout(app.PlotTab7, [1, 1], 'Padding', [8 8 8 8]);
            app.ResultsTable = uitable(tab7Grid);
            app.ResultsTable.Layout.Row = 1;
            app.ResultsTable.Layout.Column = 1;
            app.ResultsTable.ColumnName = {'Parameter', 'Value'};
            app.ResultsTable.ColumnWidth = {220, 'auto'};
            app.ResultsTable.Data = cell(0, 2);
            
            % -----------------------------------------------------------
            % BOTTOM PANEL LOG
            % -----------------------------------------------------------
            app.BottomPanel = uipanel(app.MainGridLayout);
            app.BottomPanel.Title = 'Simulation Console Log';
            app.BottomPanel.FontWeight = 'bold';
            app.BottomPanel.Layout.Row = 2;
            app.BottomPanel.Layout.Column = 2;
            
            bottomGrid = uigridlayout(app.BottomPanel, [1, 1], 'Padding', [6 6 6 6]);
            
            app.ConsoleTextArea = uitextarea(bottomGrid);
            app.ConsoleTextArea.Layout.Row = 1;
            app.ConsoleTextArea.Layout.Column = 1;
            app.ConsoleTextArea.Editable = 'off';
            app.ConsoleTextArea.BackgroundColor = [0.08 0.08 0.10]; 
            app.ConsoleTextArea.FontColor = [0.4 1.0 0.4];          
            app.ConsoleTextArea.FontName = 'Monospaced';
            app.ConsoleTextArea.FontSize = 11;
            app.ConsoleTextArea.Value = {'System ready. Select a model file and click ''Run Sim''.'};
            
            app.UIFigure.Visible = 'on';
        end
    end
    
    methods (Access = public)
        function app = app
            createComponents(app)
            registerApp(app, app.UIFigure)
            if nargout == 0
                clear app
            end
        end
        function delete(app)
            delete(app.UIFigure)
        end
    end
end