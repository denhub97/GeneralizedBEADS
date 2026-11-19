% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
classdef GeneralizedBEADSApp < matlab.apps.AppBase

    %!!! Graphical user interface for the BEADS visualization of mixed
    %    quantum states. The app is a front end for the two user functions
    %    visualize (mode "Single State") and stacked_plot (mode
    %    "Parameterized State") in the directory code. All code files are
    %    used as they are, the app only collects the inputs, calls the
    %    functions, and displays the figure files which were created.
    %
    %    Start the app from MATLAB with
    %        GeneralizedBEADSApp
    %    (the directory app must be the current folder or on the path, the
    %    directory code is added to the path automatically).
    %
    %    Inputs of the app:
    %        State: a density operator (matrix or text expression in Dirac
    %               notation, see parse) for "Single State", or a text
    %               expression / function handle (@(z) ...) of the
    %               parameter for "Parameterized State"
    %        Output: folder, file name, and format of the figures. The
    %                function appends the bead type etc. to the file name
    %        Advanced options: all optional inputs of the function of the
    %                          current mode. (default) or an empty field
    %                          means that the default value of the function
    %                          is used
    %
    %    The figures which are created by the function are loaded into the
    %    figure viewer after the computation, with buttons to navigate
    %    between the figures and to zoom. PDF files cannot be previewed.
    %
    %    Requires MATLAB R2020b or newer (as the functions in code)

    %% Components
    properties (Access = public)
        UIFigure              matlab.ui.Figure
        MainGrid              matlab.ui.container.GridLayout

        % Left panel: mode, state, output
        LeftPanel             matlab.ui.container.Panel
        LeftGrid              matlab.ui.container.GridLayout
        ModeDropDown          matlab.ui.control.DropDown
        StateLabel            matlab.ui.control.Label
        StateArea             matlab.ui.control.TextArea
        CheckButton           matlab.ui.control.Button
        CheckLabel            matlab.ui.control.Label
        FolderField           matlab.ui.control.EditField
        BrowseButton          matlab.ui.control.Button
        NameField             matlab.ui.control.EditField
        FormatDropDown        matlab.ui.control.DropDown
        PathLabel             matlab.ui.control.Label
        RunButton             matlab.ui.control.Button
        AdvancedButton        matlab.ui.control.StateButton
        StatusLabel           matlab.ui.control.Label

        % Center panel: figure viewer
        CenterPanel           matlab.ui.container.Panel
        PrevButton            matlab.ui.control.Button
        FileDropDown          matlab.ui.control.DropDown
        NextButton            matlab.ui.control.Button
        ZoomOutButton         matlab.ui.control.Button
        ZoomLabel             matlab.ui.control.Label
        ZoomInButton          matlab.ui.control.Button
        FitButton             matlab.ui.control.Button
        OpenButton            matlab.ui.control.Button
        ViewerAxes            matlab.ui.control.UIAxes
        InfoLabel             matlab.ui.control.Label

        % Right panel: advanced options
        RightPanel            matlab.ui.container.Panel
        VisGrid               matlab.ui.container.GridLayout
        StackGrid             matlab.ui.container.GridLayout
    end

    %% Data
    properties (Access = private)
        % Option definitions and the corresponding controls of both modes
        Specs = struct('single', [], 'param', [])
        Controls = struct('single', {{}}, 'param', {{}})

        % Check boxes of the bead types of both modes and the label of the
        % bead types option of the mode "Single State"
        BeadChecks = struct('single', {{}}, 'param', {{}})
        BeadLabel

        % State texts of both modes (kept when switching the mode)
        StateTexts = struct('single', '', 'param', '')

        % Figures shown in the viewer
        Files = strings(0,1)
        CurrentIndex = 0
        ImageWidth = NaN
        ImageHeight = NaN

        % Directory with the code
        CodeDir = ''
    end

    %% Component creation
    methods (Access = private)

        function createComponents(app)

            % Main window with three columns: input, viewer, options
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [60 60 1400 860];
            app.UIFigure.Name = 'Generalized BEADS';

            app.MainGrid = uigridlayout(app.UIFigure, [1 3]);
            app.MainGrid.ColumnWidth = {340, '1x', 380};
            app.MainGrid.RowHeight = {'1x'};
            app.MainGrid.ColumnSpacing = 8;
            app.MainGrid.Padding = [8 8 8 8];

            createLeftPanel(app);
            createCenterPanel(app);
            createRightPanel(app);

            app.UIFigure.Visible = 'on';
        end

        function createLeftPanel(app)

            app.LeftPanel = uipanel(app.MainGrid, 'Title', 'Input and output');
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            app.LeftGrid = uigridlayout(app.LeftPanel, [12 1]);
            app.LeftGrid.RowHeight = {26, 58, 110, 26, 48, 22, 26, 26, 44, 36, 28, '1x'};
            app.LeftGrid.RowSpacing = 6;

            % Mode selection
            modeGrid = uigridlayout(app.LeftGrid, [1 2]);
            modeGrid.Layout.Row = 1;
            modeGrid.ColumnWidth = {50, '1x'};
            modeGrid.RowHeight = {'1x'};
            modeGrid.Padding = [0 0 0 0];
            uilabel(modeGrid, 'Text', 'Mode', 'FontWeight', 'bold');
            app.ModeDropDown = uidropdown(modeGrid);
            app.ModeDropDown.Items = {'Single State', 'Parameterized State'};
            app.ModeDropDown.ItemsData = {'single', 'param'};
            app.ModeDropDown.Value = 'single';
            app.ModeDropDown.Tooltip = 'Single State: visualize, Parameterized State: stacked_plot';
            app.ModeDropDown.ValueChangedFcn = createCallbackFcn(app, @ModeChanged, true);

            % State input
            app.StateLabel = uilabel(app.LeftGrid, 'FontWeight', 'bold', 'WordWrap', 'on');
            app.StateLabel.Layout.Row = 2;
            app.StateLabel.VerticalAlignment = 'top';
            app.StateArea = uitextarea(app.LeftGrid);
            app.StateArea.Layout.Row = 3;
            app.StateArea.FontName = 'Courier New';
            app.StateArea.ValueChangedFcn = createCallbackFcn(app, @StateChanged, true);

            app.CheckButton = uibutton(app.LeftGrid, 'push', 'Text', 'Check input');
            app.CheckButton.Layout.Row = 4;
            app.CheckButton.Tooltip = 'Parse the input and check whether it is a valid density operator';
            app.CheckButton.ButtonPushedFcn = createCallbackFcn(app, @CheckButtonPushed, true);

            app.CheckLabel = uilabel(app.LeftGrid, 'Text', '', 'WordWrap', 'on');
            app.CheckLabel.Layout.Row = 5;
            app.CheckLabel.VerticalAlignment = 'top';

            % Output location
            outLabel = uilabel(app.LeftGrid, 'Text', 'Output (folder, file name, format)', 'FontWeight', 'bold');
            outLabel.Layout.Row = 6;

            folderGrid = uigridlayout(app.LeftGrid, [1 2]);
            folderGrid.Layout.Row = 7;
            folderGrid.ColumnWidth = {'1x', 80};
            folderGrid.RowHeight = {'1x'};
            folderGrid.Padding = [0 0 0 0];
            app.FolderField = uieditfield(folderGrid, 'text');
            app.FolderField.Value = pwd;
            app.FolderField.Tooltip = 'Folder in which the figures are saved';
            app.FolderField.ValueChangedFcn = createCallbackFcn(app, @OutputChanged, true);
            app.BrowseButton = uibutton(folderGrid, 'push', 'Text', 'Browse...');
            app.BrowseButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseButtonPushed, true);

            nameGrid = uigridlayout(app.LeftGrid, [1 2]);
            nameGrid.Layout.Row = 8;
            nameGrid.ColumnWidth = {'1x', 80};
            nameGrid.RowHeight = {'1x'};
            nameGrid.Padding = [0 0 0 0];
            app.NameField = uieditfield(nameGrid, 'text');
            app.NameField.Value = 'beads';
            app.NameField.Tooltip = 'File name without extension, the bead type etc. is appended by the function';
            app.NameField.ValueChangedFcn = createCallbackFcn(app, @OutputChanged, true);
            app.FormatDropDown = uidropdown(nameGrid);
            app.FormatDropDown.Items = {'png', 'jpg', 'pdf', 'tiff'};
            app.FormatDropDown.Value = 'png';
            app.FormatDropDown.Tooltip = 'Figure format (PDF files cannot be previewed)';
            app.FormatDropDown.ValueChangedFcn = createCallbackFcn(app, @OutputChanged, true);

            app.PathLabel = uilabel(app.LeftGrid, 'Text', '', 'WordWrap', 'on');
            app.PathLabel.Layout.Row = 9;
            app.PathLabel.VerticalAlignment = 'top';
            app.PathLabel.FontColor = [0.4 0.4 0.4];

            % Run and options toggle
            app.RunButton = uibutton(app.LeftGrid, 'push', 'Text', 'Generate figures');
            app.RunButton.Layout.Row = 10;
            app.RunButton.FontWeight = 'bold';
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);

            app.AdvancedButton = uibutton(app.LeftGrid, 'state', 'Text', 'Advanced options');
            app.AdvancedButton.Layout.Row = 11;
            app.AdvancedButton.Value = true;
            app.AdvancedButton.ValueChangedFcn = createCallbackFcn(app, @AdvancedButtonChanged, true);

            app.StatusLabel = uilabel(app.LeftGrid, 'Text', 'Ready.', 'WordWrap', 'on');
            app.StatusLabel.Layout.Row = 12;
            app.StatusLabel.VerticalAlignment = 'top';
        end

        function createCenterPanel(app)

            app.CenterPanel = uipanel(app.MainGrid, 'Title', 'Figure viewer');
            app.CenterPanel.Layout.Row = 1;
            app.CenterPanel.Layout.Column = 2;

            centerGrid = uigridlayout(app.CenterPanel, [3 1]);
            centerGrid.RowHeight = {28, '1x', 22};
            centerGrid.RowSpacing = 6;

            % Navigation and zoom bar
            navGrid = uigridlayout(centerGrid, [1 8]);
            navGrid.Layout.Row = 1;
            navGrid.ColumnWidth = {36, '1x', 36, 32, 56, 32, 44, 100};
            navGrid.RowHeight = {'1x'};
            navGrid.Padding = [0 0 0 0];
            navGrid.ColumnSpacing = 4;

            app.PrevButton = uibutton(navGrid, 'push', 'Text', '<');
            app.PrevButton.Tooltip = 'Previous figure';
            app.PrevButton.ButtonPushedFcn = createCallbackFcn(app, @PrevButtonPushed, true);

            app.FileDropDown = uidropdown(navGrid);
            app.FileDropDown.Items = {'(no figures)'};
            app.FileDropDown.ItemsData = {0};
            app.FileDropDown.Value = 0;
            app.FileDropDown.ValueChangedFcn = createCallbackFcn(app, @FileSelected, true);

            app.NextButton = uibutton(navGrid, 'push', 'Text', '>');
            app.NextButton.Tooltip = 'Next figure';
            app.NextButton.ButtonPushedFcn = createCallbackFcn(app, @NextButtonPushed, true);

            app.ZoomOutButton = uibutton(navGrid, 'push', 'Text', '-');
            app.ZoomOutButton.Tooltip = 'Zoom out';
            app.ZoomOutButton.ButtonPushedFcn = createCallbackFcn(app, @ZoomOutButtonPushed, true);

            app.ZoomLabel = uilabel(navGrid, 'Text', '--', 'HorizontalAlignment', 'center');

            app.ZoomInButton = uibutton(navGrid, 'push', 'Text', '+');
            app.ZoomInButton.Tooltip = 'Zoom in';
            app.ZoomInButton.ButtonPushedFcn = createCallbackFcn(app, @ZoomInButtonPushed, true);

            app.FitButton = uibutton(navGrid, 'push', 'Text', 'Fit');
            app.FitButton.Tooltip = 'Show the complete figure';
            app.FitButton.ButtonPushedFcn = createCallbackFcn(app, @FitButtonPushed, true);

            app.OpenButton = uibutton(navGrid, 'push', 'Text', 'Open externally');
            app.OpenButton.Tooltip = 'Open the figure file with the default application';
            app.OpenButton.ButtonPushedFcn = createCallbackFcn(app, @OpenButtonPushed, true);

            % Axes for the figure image (mouse wheel zoom and drag panning
            % are available in addition to the buttons)
            app.ViewerAxes = uiaxes(centerGrid);
            app.ViewerAxes.Layout.Row = 2;
            app.ViewerAxes.Visible = 'off';
            try
                app.ViewerAxes.Interactions = [zoomInteraction panInteraction];
                app.ViewerAxes.Toolbar.Visible = 'off';
            catch
                % Older releases do not support these properties
            end

            app.InfoLabel = uilabel(centerGrid, 'Text', '');
            app.InfoLabel.Layout.Row = 3;
            app.InfoLabel.FontColor = [0.4 0.4 0.4];
        end

        function createRightPanel(app)

            app.RightPanel = uipanel(app.MainGrid, 'Title', 'Advanced options (default value)');
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 3;

            % One scrollable grid with the options of each mode, only the
            % grid of the current mode is visible
            app.Specs.single = visualizeSpecs();
            app.Specs.param = stackedPlotSpecs();
            app.VisGrid = createOptionGrid(app, app.RightPanel, 'single');
            app.StackGrid = createOptionGrid(app, app.RightPanel, 'param');
            app.StackGrid.Visible = 'off';
        end

        function grid = createOptionGrid(app, parent, mode)

            % Build the option rows (label + control) from the option
            % definitions of the mode
            specs = app.Specs.(mode);
            n = numel(specs);
            grid = uigridlayout(parent, [n 2]);
            grid.ColumnWidth = {200, '1x'};
            grid.Scrollable = 'on';
            grid.RowSpacing = 4;
            rowHeights = repmat({24}, 1, n);
            controls = cell(1, n);

            for k = 1:n
                s = specs(k);

                lbl = uilabel(grid, 'Text', s.Label, 'Tooltip', s.Tip, 'WordWrap', 'on');
                lbl.Layout.Row = k;
                lbl.Layout.Column = 1;

                % Long labels need two lines
                if numel(s.Label) > 30
                    rowHeights{k} = 36;
                end

                switch s.Type
                    case 'bool'
                        c = uidropdown(grid, 'Items', {'(default)', 'true', 'false'}, 'Value', '(default)');
                    case 'enum'
                        c = uidropdown(grid, 'Items', [{'(default)'}, s.Items], 'Value', '(default)');
                    case 'structure'
                        % Structure of the optimization, the items are
                        % mapped to -1, 0, 1 in collectOptions
                        c = uidropdown(grid, 'Items', {'entangling', 'none (default)', 'separable'}, 'Value', 'none (default)');
                    case 'multi'
                        % One check box per item in two columns
                        nRows = ceil(numel(s.Items) / 2);
                        c = uigridlayout(grid, [nRows 2]);
                        c.Padding = [0 0 0 0];
                        c.RowSpacing = 2;
                        c.ColumnSpacing = 4;
                        c.RowHeight = repmat({22}, 1, nRows);
                        checks = cell(1, numel(s.Items));
                        for m = 1:numel(s.Items)
                            checks{m} = uicheckbox(c, 'Text', s.Items{m}, 'Value', false, 'Tooltip', s.Tip);
                        end
                        app.BeadChecks.(mode) = checks;
                        rowHeights{k} = 24 * nRows + 4;
                        if strcmp(mode, 'single')
                            app.BeadLabel = lbl;
                        end
                    otherwise
                        c = uieditfield(grid, 'text', 'Value', '');
                end
                c.Layout.Row = k;
                c.Layout.Column = 2;
                if ~strcmp(s.Type, 'multi')
                    c.Tooltip = s.Tip;
                end

                % The overview mode depends on the SingleView option
                if strcmp(mode, 'single') && strcmp(s.Name, 'SingleView')
                    c.ValueChangedFcn = createCallbackFcn(app, @SingleViewChanged, true);
                end
                controls{k} = c;
            end

            grid.RowHeight = rowHeights;
            app.Controls.(mode) = controls;
        end
    end

    %% Helper methods
    methods (Access = private)

        function bringToFront(app)

            % Bring the app window back to the front after other figures
            % or dialogs were shown
            try
                figure(app.UIFigure);
            catch
                % Not supported in this release, nothing to do
            end
        end

        function setStatus(app, msg, color)

            % Show a message in the status label
            if nargin < 3
                color = [0.2 0.2 0.2];
            end
            app.StatusLabel.Text = msg;
            app.StatusLabel.FontColor = color;
            drawnow;
        end

        function setCheck(app, msg, color)

            % Show the result of the input check
            app.CheckLabel.Text = msg;
            app.CheckLabel.FontColor = color;
        end

        function updateOptionAvailability(app)

            % Block options which have no effect for the current input. In
            % the mode "Single State": for a two-qubit state, visualize
            % creates one overview figure of all bead types if SingleView
            % is false or default, then the BeadTypes option has no effect,
            % and the OptTol option only applies to three qubits. The
            % Detailed option only applies to the overview figure and is
            % blocked if SingleView is true. Blocked options are not
            % passed to the function
            checks = app.BeadChecks.single;
            if isempty(checks)
                return;
            end

            % Number of qubits of the input state (0 if it is invalid)
            nQubits = 0;
            if strcmp(app.ModeDropDown.Value, 'single')
                try
                    text = getStateText(app);
                    if text ~= ""
                        nQubits = log2(size(parse(text), 1));
                    end
                catch
                    % Invalid input, the options stay available
                end
            end

            names = {app.Specs.single.Name};
            svValue = app.Controls.single{find(strcmp(names, 'SingleView'), 1)}.Value;
            overview = nQubits == 2 && any(strcmp(svValue, {'(default)', 'false'}));

            % Bead types
            if overview
                state = 'off';
                app.BeadLabel.Text = 'Bead types (overview: all)';
            else
                state = 'on';
                app.BeadLabel.Text = 'Bead types (all)';
            end
            for m = 1:numel(checks)
                checks{m}.Enable = state;
            end

            % Optimization tolerance
            if nQubits == 2
                state = 'off';
            else
                state = 'on';
            end
            app.Controls.single{find(strcmp(names, 'OptTol'), 1)}.Enable = state;

            % Detailed set of beads, only applies to the overview figure
            % (SingleView false or default)
            if strcmp(svValue, 'true')
                state = 'off';
            else
                state = 'on';
            end
            app.Controls.single{find(strcmp(names, 'Detailed'), 1)}.Enable = state;
        end

        function updateViewerButtons(app)

            % Enable the buttons of the viewer only if they can be used
            nFiles = numel(app.Files);
            hasImage = ~isnan(app.ImageWidth);
            onoff = {'off', 'on'};
            app.PrevButton.Enable = onoff{(nFiles > 1) + 1};
            app.NextButton.Enable = onoff{(nFiles > 1) + 1};
            app.FileDropDown.Enable = onoff{(nFiles > 0) + 1};
            app.ZoomOutButton.Enable = onoff{hasImage + 1};
            app.ZoomInButton.Enable = onoff{hasImage + 1};
            app.FitButton.Enable = onoff{hasImage + 1};
            app.OpenButton.Enable = onoff{(nFiles > 0) + 1};
        end

        function text = getStateText(app)

            % Get the state input as a single line of text
            text = strtrim(strjoin(string(app.StateArea.Value), " "));
        end

        function [path, folder, base, ext] = getOutputPath(app)

            % Assemble and validate the output path from folder, file name,
            % and format
            folder = strtrim(app.FolderField.Value);
            assert(isfolder(folder), "The output folder does not exist.");
            ext = app.FormatDropDown.Value;

            % A typed extension in the file name is removed
            base = regexprep(strtrim(app.NameField.Value), '\.(png|jpg|pdf|tiff)$', '', 'ignorecase');
            assert(~isempty(base) && isempty(regexp(base, '[\\/:*?"<>|]', 'once')), ...
                   "The file name is empty or contains forbidden characters (slash, backslash, colon, asterisk, question mark, quotes, angle brackets, bar).");

            path = string(fullfile(folder, [base '.' ext]));
        end

        function args = collectOptions(app, mode)

            % Collect all options which are not left at their default as
            % name-value pairs
            specs = app.Specs.(mode);
            controls = app.Controls.(mode);
            args = {};

            for k = 1:numel(specs)
                s = specs(k);
                c = controls{k};

                % Blocked options are not passed
                if ~strcmp(s.Type, 'multi') && strcmp(c.Enable, 'off')
                    continue;
                end

                switch s.Type
                    case 'bool'
                        if ~strcmp(c.Value, '(default)')
                            args = [args, {s.Name, strcmp(c.Value, 'true')}]; %#ok<AGROW>
                        end
                    case 'enum'
                        if ~strcmp(c.Value, '(default)')
                            args = [args, {s.Name, string(c.Value)}]; %#ok<AGROW>
                        end
                    case 'structure'
                        % entangling = -1, none = 0 (default), separable = 1
                        if ~strcmp(c.Value, 'none (default)')
                            args = [args, {s.Name, double(strcmp(c.Value, 'separable')) - double(strcmp(c.Value, 'entangling'))}]; %#ok<AGROW>
                        end
                    case 'multi'
                        % Checked bead types, nothing checked = default. A
                        % blocked option (see updateOptionAvailability) is not
                        % passed
                        checks = app.BeadChecks.(mode);
                        selected = {};
                        blocked = false;
                        for m = 1:numel(checks)
                            blocked = blocked || strcmp(checks{m}.Enable, 'off');
                            if checks{m}.Value
                                selected{end+1} = checks{m}.Text; %#ok<AGROW>
                            end
                        end
                        if ~blocked && ~isempty(selected)
                            args = [args, {s.Name, string(selected)}]; %#ok<AGROW>
                        end
                    case 'text'
                        txt = strtrim(c.Value);
                        if ~isempty(txt)
                            args = [args, {s.Name, string(txt)}]; %#ok<AGROW>
                        end
                    case 'num'
                        txt = strtrim(c.Value);
                        if ~isempty(txt)
                            args = [args, {s.Name, parseNumeric(txt, s.Name)}]; %#ok<AGROW>
                        end
                end
            end
        end

        function handle = buildParamHandle(~, text, name)

            % Function handle of the parametrized state from a text, either
            % an anonymous function (starting with @) or an expression
            if startsWith(text, "@")
                handle = str2func(char(text));
            else
                handle = parse(text, "Parameter", name);
            end
        end

        function files = listOutputFiles(~, folder, base, ext)

            % Map of the files which can be created for the base name with
            % their modification times
            files = containers.Map('KeyType', 'char', 'ValueType', 'double');
            d = dir(fullfile(folder, [base '*.' ext]));
            d = d(~[d.isdir]);
            for k = 1:numel(d)
                files(d(k).name) = d(k).datenum;
            end
        end

        function loadFiles(app, files)

            % Load a list of figure files into the viewer
            app.Files = files(:);
            app.FileDropDown.ItemsData = {};
            if isempty(files)
                app.FileDropDown.Items = {'(no figures)'};
                app.FileDropDown.ItemsData = {0};
                app.FileDropDown.Value = 0;
                showEmpty(app, 'No figure files were found.');
                return;
            end

            names = cell(1, numel(files));
            for k = 1:numel(files)
                [~, n, e] = fileparts(files(k));
                names{k} = sprintf('%d/%d  %s%s', k, numel(files), n, e);
            end
            app.FileDropDown.Items = names;
            app.FileDropDown.ItemsData = num2cell(1:numel(files));
            showFile(app, 1);
        end

        function showEmpty(app, msg)

            % Show a message instead of a figure
            ax = app.ViewerAxes;
            cla(ax);
            axis(ax, 'normal');
            ax.XLim = [0 1];
            ax.YLim = [0 1];
            ax.YDir = 'normal';
            ax.Visible = 'off';
            text(ax, 0.5, 0.5, msg, 'HorizontalAlignment', 'center', 'Interpreter', 'none');
            app.ImageWidth = NaN;
            app.ImageHeight = NaN;
            app.ZoomLabel.Text = '--';
            app.InfoLabel.Text = '';
            updateViewerButtons(app);
        end

        function showFile(app, idx)

            % Show the figure with the given index in the viewer
            app.CurrentIndex = idx;
            app.FileDropDown.Value = idx;
            file = char(app.Files(idx));
            [~, ~, ext] = fileparts(file);
            app.InfoLabel.Text = file;

            % PDF files cannot be displayed
            if strcmpi(ext, '.pdf')
                showEmpty(app, "PDF files cannot be previewed." + newline + "Use 'Open externally'.");
                app.InfoLabel.Text = file;
                return;
            end

            try
                [img, map] = imread(file);
                if ~isempty(map)
                    img = ind2rgb(img, map);
                end
                if size(img, 3) == 1
                    img = repmat(img, [1 1 3]);
                elseif size(img, 3) > 3
                    img = img(:, :, 1:3);
                end
                if isa(img, 'uint16')
                    img = uint8(double(img) / 257);
                end
            catch ME
                showEmpty(app, "The figure could not be loaded:" + newline + string(ME.message));
                app.InfoLabel.Text = file;
                return;
            end

            ax = app.ViewerAxes;
            cla(ax);
            image(ax, img);
            axis(ax, 'image');
            ax.YDir = 'reverse';
            ax.Visible = 'off';
            app.ImageWidth = size(img, 2);
            app.ImageHeight = size(img, 1);
            fitView(app);
            updateViewerButtons(app);
        end

        function fitView(app)

            % Show the complete image
            if isnan(app.ImageWidth)
                return;
            end
            app.ViewerAxes.XLim = [0.5, app.ImageWidth + 0.5];
            app.ViewerAxes.YLim = [0.5, app.ImageHeight + 0.5];
            updateZoomLabel(app);
        end

        function zoomBy(app, factor)

            % Zoom about the center of the current view (limited to the
            % image and a maximum zoom of 20)
            if isnan(app.ImageWidth)
                return;
            end
            ax = app.ViewerAxes;
            cx = mean(ax.XLim);
            cy = mean(ax.YLim);

            newW = diff(ax.XLim) / factor;
            newW = min(max(newW, app.ImageWidth / 20), app.ImageWidth);
            newH = newW * app.ImageHeight / app.ImageWidth;

            % Keep the view inside the image
            x0 = min(max(cx - newW/2, 0.5), app.ImageWidth + 0.5 - newW);
            y0 = min(max(cy - newH/2, 0.5), app.ImageHeight + 0.5 - newH);
            ax.XLim = [x0, x0 + newW];
            ax.YLim = [y0, y0 + newH];
            updateZoomLabel(app);
        end

        function updateZoomLabel(app)

            % Show the current zoom level
            if isnan(app.ImageWidth)
                app.ZoomLabel.Text = '--';
            else
                app.ZoomLabel.Text = sprintf('%d%%', round(100 * app.ImageWidth / diff(app.ViewerAxes.XLim)));
            end
        end
    end

    %% Callbacks
    methods (Access = private)

        function startupFcn(app)

            % Add the code directory to the path and set the default texts
            appDir = fileparts(mfilename('fullpath'));
            app.CodeDir = fullfile(fileparts(appDir), 'Code');
            if isfolder(app.CodeDir)
                addpath(genpath(app.CodeDir));
            else
                uialert(app.UIFigure, "The directory code was not found next to the directory app: " + string(app.CodeDir), ...
                        'Code directory missing');
            end

            app.StateTexts.single = '0.7*dens(|Phi+>) + 0.3*eye(4)/4';
            app.StateTexts.param = 'z*dens(|Phi+>) + (1-z)*eye(4)/4';
            app.StateArea.Value = app.StateTexts.single;
            ModeChanged(app, struct('Value', 'single', 'PreviousValue', 'single'));
            OutputChanged(app);
            showEmpty(app, 'No figure generated yet.');
            try
                movegui(app.UIFigure, 'center');
            catch
                % Keep the default position
            end
        end

        function ModeChanged(app, event)

            % Store the state text of the old mode and switch the options
            % and the state input to the new mode
            if ~strcmp(event.PreviousValue, event.Value)
                app.StateTexts.(event.PreviousValue) = app.StateArea.Value;
            end
            mode = event.Value;
            app.StateArea.Value = app.StateTexts.(mode);

            if strcmp(mode, 'single')
                app.StateLabel.Text = 'Density operator (matrix or Dirac expression, e.g. |Phi+><Phi+|):';
                app.VisGrid.Visible = 'on';
                app.StackGrid.Visible = 'off';
            else
                app.StateLabel.Text = 'Parametrized state (expression of the parameter, or @(z) ...). The parameter name is set in the options (default z):';
                app.VisGrid.Visible = 'off';
                app.StackGrid.Visible = 'on';
            end
            setCheck(app, '', [0.2 0.2 0.2]);
            updateOptionAvailability(app);
        end

        function StateChanged(app, ~)

            % The state input changed, the check result is outdated
            setCheck(app, '', [0.2 0.2 0.2]);
            updateOptionAvailability(app);
        end

        function SingleViewChanged(app, ~)

            updateOptionAvailability(app);
        end

        function OutputChanged(app, ~)

            % Show the resulting base path
            try
                path = getOutputPath(app);
                app.PathLabel.Text = "Saved as: " + path;
                app.PathLabel.FontColor = [0.4 0.4 0.4];
            catch ME
                app.PathLabel.Text = string(ME.message);
                app.PathLabel.FontColor = [0.75 0 0];
            end
        end

        function BrowseButtonPushed(app, ~)

            % Select the output folder
            start = app.FolderField.Value;
            if ~isfolder(start)
                start = pwd;
            end
            folder = uigetdir(start, 'Select the output folder');
            bringToFront(app);
            if ischar(folder) || isstring(folder)
                app.FolderField.Value = char(folder);
                OutputChanged(app);
            end
        end

        function AdvancedButtonChanged(app, ~)

            % Show or hide the advanced options
            % Only the column width is changed, switching the visibility of
            % the panel with all its controls is slow
            if app.AdvancedButton.Value
                app.MainGrid.ColumnWidth = {340, '1x', 380};
            else
                app.MainGrid.ColumnWidth = {340, '1x', 0};
            end
        end

        function CheckButtonPushed(app, ~)

            % Parse the input and report whether it is a valid state
            mode = app.ModeDropDown.Value;
            updateOptionAvailability(app);
            text = getStateText(app);
            if text == ""
                setCheck(app, 'The input is empty.', [0.75 0 0]);
                return;
            end

            try
                if strcmp(mode, 'single')
                    rho = parse(text);
                    setCheck(app, sprintf('Valid %d-qubit density operator (rank %d, purity %.4f).', ...
                                          log2(size(rho,1)), rank(rho), real(trace(rho*rho))), [0 0.5 0]);
                else
                    args = collectOptions(app, 'param');
                    name = char(strtrim(erase(string(argValue(args, 'Parameter', "z")), "$")));
                    values = argValue(args, 'ParamValues', 0:0.2:1);
                    handle = buildParamHandle(app, text, name);

                    % Evaluate the state for the first and the last value
                    for v = [values(1), values(end)]
                        rho = handle(v);
                        assert(isnumeric(rho) && isequal(size(rho), [4 4]), ...
                               "The state must be a 4x4 density operator (two qubits), but is %d x %d for %s = %g.", ...
                               size(rho,1), size(rho,2), name, v);
                    end
                    setCheck(app, sprintf('Valid two-qubit state for %s = %g and %s = %g.', name, values(1), name, values(end)), [0 0.5 0]);
                end
            catch ME
                setCheck(app, ME.message, [0.75 0 0]);
            end
        end

        function RunButtonPushed(app, ~)

            % Collect and validate all inputs
            mode = app.ModeDropDown.Value;
            updateOptionAvailability(app);
            try
                [path, folder, base, ext] = getOutputPath(app);
                text = getStateText(app);
                assert(text ~= "", "The state input is empty.");
                args = collectOptions(app, mode);
                if strcmp(mode, 'param') && startsWith(text, "@")
                    state = str2func(char(text));
                else
                    state = text;
                end
            catch ME
                uialert(app.UIFigure, ME.message, 'Invalid input');
                return;
            end

            % Ask for confirmation if existing files may be overwritten. The
            % functions append the bead type etc. to the file name, which
            % files are written depends on the options
            before = listOutputFiles(app, folder, base, ext);
            existing = string(keys(before));
            if ~isempty(existing)
                pattern = ['^' regexptranslate('escape', base) '(_.*)?\.' ext '$'];
                existing = existing(~cellfun(@isempty, regexp(cellstr(existing), pattern, 'once', 'ignorecase')));
            end
            if ~isempty(existing)
                shown = existing(1:min(numel(existing), 8));
                msg = "The following existing files in the output folder may be overwritten:" + newline + newline + strjoin(shown, newline);
                if numel(existing) > numel(shown)
                    msg = msg + newline + "... and " + (numel(existing) - numel(shown)) + " more";
                end
                msg = msg + newline + newline + "Continue?";
                choice = uiconfirm(app.UIFigure, msg, 'Overwrite existing files?', ...
                                   'Options', {'Overwrite', 'Cancel'}, 'DefaultOption', 2, ...
                                   'CancelOption', 2, 'Icon', 'warning');
                if strcmp(choice, 'Cancel')
                    return;
                end
            end

            % The app never uses the prompt of the code for the measurement
            % direction, the z-axis is the default
            if ~any(strcmp(args(1:2:end), 'MeasurementDirection'))
                args = [args, {'MeasurementDirection', [0 0 1]}];
            end

            % Remember the existing figure windows
            figsBefore = findall(0, 'Type', 'figure');

            % Run the function of the mode
            app.RunButton.Enable = 'off';
            setStatus(app, 'Computing... (progress is shown in the command window)');
            dlg = uiprogressdlg(app.UIFigure, 'Title', 'Generating figures', ...
                                'Message', 'Computing the beads, this can take a while.', 'Indeterminate', 'on');
            success = true;
            try
                if strcmp(mode, 'single')
                    visualize(state, path, args{:});
                else
                    stacked_plot(state, path, args{:});
                end
            catch ME
                success = false;
                msg = ME.message;
            end

            % Close the progress dialog and leftover figure windows of an
            % aborted run
            close(dlg);
            if ~success
                figsNow = findall(0, 'Type', 'figure');
                delete(setdiff(figsNow, figsBefore));
            end
            bringToFront(app);
            app.RunButton.Enable = 'on';

            if ~success
                setStatus(app, 'The computation failed.', [0.75 0 0]);
                uialert(app.UIFigure, msg, 'Error while generating figures');
                return;
            end

            % Find the files which were created or changed by the run
            after = listOutputFiles(app, folder, base, ext);
            names = string(keys(after));
            isNew = false(size(names));
            for k = 1:numel(names)
                isNew(k) = ~isKey(before, char(names(k))) || after(char(names(k))) > before(char(names(k)));
            end
            names = sort(names(isNew));
            if isempty(names)
                files = strings(0,1);
            else
                files = string(fullfile(folder, names));
            end

            loadFiles(app, files);
            if isempty(files)
                setStatus(app, 'The computation finished, but no new figure files were found.', [0.75 0.4 0]);
            else
                setStatus(app, sprintf('Done: %d figure(s) saved in %s', numel(files), folder), [0 0.5 0]);
            end
        end

        function PrevButtonPushed(app, ~)

            % Previous figure
            if app.CurrentIndex > 1
                showFile(app, app.CurrentIndex - 1);
            end
        end

        function NextButtonPushed(app, ~)

            % Next figure
            if app.CurrentIndex >= 1 && app.CurrentIndex < numel(app.Files)
                showFile(app, app.CurrentIndex + 1);
            end
        end

        function FileSelected(app, event)

            % Figure selected in the drop-down list
            if event.Value >= 1
                showFile(app, event.Value);
            end
        end

        function ZoomInButtonPushed(app, ~)

            zoomBy(app, 1.25);
        end

        function ZoomOutButtonPushed(app, ~)

            zoomBy(app, 1/1.25);
        end

        function FitButtonPushed(app, ~)

            fitView(app);
        end

        function OpenButtonPushed(app, ~)

            % Open the current figure file with the default application
            if app.CurrentIndex < 1
                return;
            end
            file = char(app.Files(app.CurrentIndex));
            if ispc
                winopen(file);
            elseif ismac
                system(['open "' file '"']);
            else
                system(['xdg-open "' file '" &']);
            end
        end
    end

    %% App creation and deletion
    methods (Access = public)

        function app = GeneralizedBEADSApp

            % Create the components, register the app, and run the startup
            createComponents(app);
            registerApp(app, app.UIFigure);
            runStartupFcn(app, @startupFcn);

            if nargout == 0
                clear app;
            end
        end

        function delete(app)

            % Delete the window when the app is deleted
            delete(app.UIFigure);
        end
    end
end

%% Local functions

function specs = visualizeSpecs()

    %!!! Definitions of the options of visualize: name, control type,
    %    default (for the label), tooltip, and items

    beadTypes = {'T', 'E', 'C', 'M', 'Sigma', 'M_Cl', 'E_Cl', 'Sigma_Cl', 'D', 'E_D', 'Sigma_D'};
    specs = [ ...
        makeSpec('BeadTypes', 'multi', 'all', 'Bead types to plot. Nothing checked = all nine bead types.', beadTypes); ...
        makeSpec('SingleView', 'bool', 'false', 'Two qubits only: false = one overview figure with all bead types (BeadTypes is ignored), true = one figure per bead type.', {}); ...
        makeSpec('PureStateVis', 'bool', 'false', 'Also plot the BEADS of the individual pure states of the decomposition.', {}); ...
        makeSpec('Detailed', 'bool', 'false', 'Two qubits only: plot the detailed set of beads in the overview figure.', {}); ...
        makeSpec('CameraAngle', 'num', 'auto', 'Camera perspective [elevation, rotation], e.g. [47 30]. Default: [47 30] (2 qubits), [54 45] (3 qubits).', {}); ...
        makeSpec('GridPoints', 'num', '200 400', 'Resolution of the spherical grid [theta, phi], e.g. [200 400].', {}); ...
        makeSpec('TargetQubits', 'num', '1 1 2 1 2', 'Target qubits of the classical beads, e.g. [1 1 2 1 2] (for two qubits only the first entry is used).', {}); ...
        makeSpec('Id', 'bool', 'false', 'Draw the identity bead.', {}); ...
        makeSpec('EntanglementArcs', 'bool', 'true', 'Draw the entanglement arcs between the beads.', {}); ...
        makeSpec('RadialPlot', 'bool', 'false', 'Radial plot: the distance from the origin encodes the absolute function value.', {}); ...
        makeSpec('AutoScale', 'bool', 'false', 'Scale the function values relative to the global maximum.', {}); ...
        makeSpec('AddRadialPlot', 'bool', 'false', 'Additionally create a radial version of each figure.', {}); ...
        makeSpec('AutoScaleAddRadial', 'bool', 'false', 'Autoscale the additional radial plots.', {}); ...
        makeSpec('PadVectors', 'num', 'auto', 'Number of zero vectors added to the pure state decomposition. Default: rank^2 + 1 - rank. 0 switches the padding off.', {}); ...
        makeSpec('OptStructure', 'structure', 'none', 'Structure which the classical beads are optimized for: entangling = minimize the overlap with the separable structure, none = no structure optimization, separable = maximize the overlap with the separable structure (the optimizations are slow).', {}); ...
        makeSpec('MeasurementDirection', 'num', 'z-axis', 'Measurement direction [x y z] of the classical beads, used if the optimal direction is arbitrary (two qubits and marginals). The app uses [0 0 1] if empty.', {}); ...
        makeSpec('OptTol', 'num', '1e-6', 'Tolerance of the optimization of the pure state decomposition (three qubits only).', {}); ...
        makeSpec('FontSize', 'num', '15', 'Font size of all labels and annotations.', {})];
end

function specs = stackedPlotSpecs()

    %!!! Definitions of the options of stacked_plot: name, control type,
    %    default (for the label), tooltip, and items

    beadTypes = {'M', 'T', 'E', 'C', 'Sigma', 'Sigma_Cl', 'Sigma_D', 'D', 'M_Cl'};
    specs = [ ...
        makeSpec('Parameter', 'text', 'z', 'Name of the parameter in the state expression (without dollar signs), plotted as $name$.', {}); ...
        makeSpec('ParamValues', 'num', '0:0.2:1', 'Parameter values of the stacked states, e.g. 0:0.25:1 or [0 0.5 1].', {}); ...
        makeSpec('BeadTypes', 'multi', 'all', 'Bead types to plot, one figure per type. Nothing checked = all nine bead types.', beadTypes); ...
        makeSpec('Id', 'bool', 'false', 'Draw the identity bead.', {}); ...
        makeSpec('EntArcs', 'bool', 'true', 'Draw the connector arcs, scaled by the entanglement of each state.', {}); ...
        makeSpec('EntMeasure', 'enum', 'tangle', 'Entanglement measure for the connectors: two-tangle or entanglement of formation.', {'tangle', 'EoF'}); ...
        makeSpec('GridPoints', 'num', '200 400', 'Resolution of the spherical grid [theta, phi], e.g. [200 400].', {}); ...
        makeSpec('OptStructure', 'structure', 'none', 'Structure which the classical beads are optimized for: entangling = minimize the overlap with the separable structure, none = no structure optimization, separable = maximize the overlap with the separable structure (the optimizations are slow).', {}); ...
        makeSpec('MeasurementDirection', 'num', 'z-axis', 'Measurement direction [x y z] of the classical beads, used if the optimal direction is arbitrary. The app uses [0 0 1] if empty.', {}); ...
        makeSpec('TargetQubit', 'num', '1', 'Measured qubit (1 or 2) for the classical beads.', {}); ...
        makeSpec('FontSize', 'num', '12', 'Font size of all labels and annotations.', {})];
end

function spec = makeSpec(name, type, defaultText, tip, items)

    %!!! Create the definition of one option. The label shows the display
    %    name of the option (see displayName) and the default value

    spec = struct('Name', name, 'Type', type, 'Label', sprintf('%s (%s)', displayName(name), defaultText), ...
                  'Tip', sprintf('%s Default: %s.', tip, defaultText), 'Items', {items});
end

function value = parseNumeric(text, name)

    %!!! Convert the text of a numerical option into numbers. Only digits,
    %    the characters of matrices, ranges and arithmetic expressions, and
    %    the words sqrt and pi are allowed, e.g., "[200 400]", "0:0.25:1",
    %    or "[1/sqrt(2) 0 1/sqrt(2)]"

    allowedOnly = regexprep(text, '\<(sqrt|pi)\>', '');
    if ~isempty(regexp(allowedOnly, '[^0-9eE\s\[\],;:.+\-*/^()]', 'once'))
        error("The option %s must be a number or a vector, e.g. [200 400], 0:0.25:1, or [1/sqrt(2) 0 1/sqrt(2)].", name);
    end
    [value, ok] = str2num(text); %#ok<ST2NM>
    if ~ok || ~isnumeric(value) || isempty(value)
        error("The option %s could not be interpreted as numbers.", name);
    end
end

function value = argValue(args, name, default)

    %!!! Get the value of a name-value pair from a cell array of
    %    arguments, or the default if it is not contained

    idx = find(strcmp(args(1:2:end), name), 1);
    if isempty(idx)
        value = default;
    else
        value = args{2*idx};
    end
end

function label = displayName(name)

    %!!! Display name of an option in the advanced options menu, the
    %    names of the options themselves are not changed

    switch name
        case 'BeadTypes',            label = 'Bead types';
        case 'SingleView',           label = 'Single view';
        case 'PureStateVis',         label = 'Visualize pure states';
        case 'Detailed',             label = 'Plot detailed set of beads';
        case 'CameraAngle',          label = 'Camera angle';
        case 'GridPoints',           label = 'Grid points';
        case 'TargetQubits',         label = 'Target qubits';
        case 'TargetQubit',          label = 'Target qubit';
        case 'Id',                   label = 'Show Id-Bead';
        case {'EntanglementArcs', 'EntArcs'}, label = 'Entanglement arcs';
        case 'EntMeasure',           label = 'Entanglement measure';
        case 'RadialPlot',           label = 'Radial plot';
        case 'AutoScale',            label = 'Autoscale beads';
        case 'AddRadialPlot',        label = 'Additional radial plot';
        case 'AutoScaleAddRadial',   label = 'Autoscale beads of radial plot';
        case 'PadVectors',           label = 'Added zero vectors for optim.';
        case 'OptStructure',         label = 'Maximize structural contribution of M_Cl-Beads';
        case 'MeasurementDirection', label = 'Default measurement direction';
        case 'OptTol',               label = 'Optimization tolerance';
        case 'FontSize',             label = 'Font size';
        case 'Parameter',            label = 'Parameter name';
        case 'ParamValues',          label = 'Parameter values';
        otherwise,                   label = name;
    end
end
