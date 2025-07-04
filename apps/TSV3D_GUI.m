%% 3D-TSV
%%The 3D Trajectory-based Stress Visualizer (3D-TSV), a visual analysis tool for the exploration 
%%of the principal stress directions in 3D solids under load.

%%This repository was created for the paper "3D-TSV: The 3D Trajectory-based Stress Visualizer" 
%%	by Junpeng Wang, Christoph Neuhauser, Jun Wu, Xifeng Gao and Rüdiger Westermann, 
%%which was submitted to the journal "Advances in Engineering Software", and also available in arXiv (2112.09202)

%% ============= Test data sets can be found in =============
%%	https://syncandshare.lrz.de/getlink/fi4W4EGjZSzMzCvxkEf9L3Aw/

classdef TSV3D_GUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        GridLayout                      matlab.ui.container.GridLayout
        LeftPanel                       matlab.ui.container.Panel
        SubmitButton                    matlab.ui.control.Button
        TabGroup                        matlab.ui.container.TabGroup
        SettingsTab                     matlab.ui.container.Tab
        SnappingCheckBox                matlab.ui.control.CheckBox
        MergingCheckBox                 matlab.ui.control.CheckBox
        MinorCheckBox                   matlab.ui.control.CheckBox
        MediumCheckBox                  matlab.ui.control.CheckBox
        MajorCheckBox                   matlab.ui.control.CheckBox
        SeedingStrategyDropDown         matlab.ui.control.DropDown
        SeedingStrategyDropDownLabel    matlab.ui.control.Label
        PermittedMaxTangentDeviEditFieldLabel  matlab.ui.control.Label
        PermittedMaxTangentDeviEditField  matlab.ui.control.NumericEditField
        IntegratingSchemeDropDownLabel  matlab.ui.control.Label
        IntegratingSchemeDropDown       matlab.ui.control.DropDown
        MergingThresholdScalingFacMajorEditFieldLabel  matlab.ui.control.Label
        MergingThresholdScalingFacMajorEditField  matlab.ui.control.NumericEditField
        MergingThresholdScalingFacMediumEditFieldLabel  matlab.ui.control.Label
        MergingThresholdScalingFacMediumEditField  matlab.ui.control.NumericEditField
        MergingThresholdScalingFacMinorEditFieldLabel  matlab.ui.control.Label
        MergingThresholdScalingFacMinorEditField  matlab.ui.control.NumericEditField
        PSLGeometryScalingFacLabel      matlab.ui.control.Label
        PSLGeometryScalingFacEditField  matlab.ui.control.NumericEditField
        PermittedMinLengthofVisiblePSLEditFieldLabel  matlab.ui.control.Label
        PermittedMinLengthofVisiblePSLEditField  matlab.ui.control.NumericEditField
        VoxelizedPSLThicknessExportEditFieldLabel  matlab.ui.control.Label
        VoxelizedPSLThicknessExportEditField  matlab.ui.control.NumericEditField
        BoundaryThicknessExportEditFieldLabel  matlab.ui.control.Label
        BoundaryThicknessExportEditField  matlab.ui.control.NumericEditField
        SeedDensityCtrlEditField_2Label  matlab.ui.control.Label
        SeedDensityCtrlEditField_2      matlab.ui.control.NumericEditField
        StepsizeScalingFacdelta_sEditFieldLabel  matlab.ui.control.Label
        StepsizeScalingFacdelta_sEditField  matlab.ui.control.NumericEditField
        InteractionsTab                 matlab.ui.container.Tab
        MajorLoDsDropDownLabel          matlab.ui.control.Label
        MajorLoDsDropDown               matlab.ui.control.DropDown
        MajorPSLGeometryDropDownLabel   matlab.ui.control.Label
        MajorPSLGeometryDropDown        matlab.ui.control.DropDown
        MajorResCtrlSliderLabel         matlab.ui.control.Label
        MajorResCtrlSlider              matlab.ui.control.Slider
        MediumResCtrlSliderLabel        matlab.ui.control.Label
        MediumResCtrlSlider             matlab.ui.control.Slider
        MediumLoDsDropDownLabel         matlab.ui.control.Label
        MediumLoDsDropDown              matlab.ui.control.DropDown
        MediumPSLGeometryDropDownLabel  matlab.ui.control.Label
        MediumPSLGeometryDropDown       matlab.ui.control.DropDown
        MinorResCtrlSliderLabel         matlab.ui.control.Label
        MinorResCtrlSlider              matlab.ui.control.Slider
        MinorLoDsDropDownLabel          matlab.ui.control.Label
        MinorLoDsDropDown               matlab.ui.control.DropDown
        MinorPSLGeometryDropDownLabel   matlab.ui.control.Label
        MinorPSLGeometryDropDown        matlab.ui.control.DropDown
        ColoringDropDownLabel           matlab.ui.control.Label
        ColoringDropDown                matlab.ui.control.DropDown
        CameraZoominoutEditFieldLabel   matlab.ui.control.Label
        CameraZoominoutEditField        matlab.ui.control.NumericEditField
        LineDensityCtrlEditField_2Label  matlab.ui.control.Label
        LineDensityCtrlEditField_2      matlab.ui.control.NumericEditField
        NumLevelsEditField_2Label       matlab.ui.control.Label
        NumLevelsEditField_2            matlab.ui.control.NumericEditField
        RightPanel                      matlab.ui.container.Panel
        UIAxes                          matlab.ui.control.UIAxes
        FilesMenu                       matlab.ui.container.Menu
        ImportStressFieldMenu           matlab.ui.container.Menu
        ExportPSLsforLineVisdatMenu     matlab.ui.container.Menu
        MediaMenu                       matlab.ui.container.Menu
        ProblemDescriptionMenu          matlab.ui.container.Menu
        MeshMenu                        matlab.ui.container.Menu
        ScalarStressFieldMenu           matlab.ui.container.Menu
        vonMisesStressMenu              matlab.ui.container.Menu
        NormalStressMenu                matlab.ui.container.Menu
        xxMenu                          matlab.ui.container.Menu
        yyMenu                          matlab.ui.container.Menu
        zzMenu                          matlab.ui.container.Menu
        ShearStressMenu                 matlab.ui.container.Menu
        yzMenu                          matlab.ui.container.Menu
        zxMenu                          matlab.ui.container.Menu
        xyMenu                          matlab.ui.container.Menu
        PrincipalStressMenu             matlab.ui.container.Menu
        MajorMenu                       matlab.ui.container.Menu
        MediumMenu                      matlab.ui.container.Menu
        MinorMenu                       matlab.ui.container.Menu
        SeedsMenu                       matlab.ui.container.Menu
        PSLsMenu                        matlab.ui.container.Menu
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: SubmitButton
        function SubmitButtonPushed(app, event)
            global runEnvironment_;
            global userInterface_;
            global integratingStepScalingFac_;
            if isempty(runEnvironment_), return; end
            userInterface_.lineDensCtrl = app.LineDensityCtrlEditField_2.Value;
            userInterface_.numLevels = app.NumLevelsEditField_2.Value;
            switch app.SeedingStrategyDropDown.Value
                case 'Volume', userInterface_.seedStrategy = 'Volume';
                case 'Surface', userInterface_.seedStrategy = 'Surface';
                case 'Loading Area', userInterface_.seedStrategy = 'LoadingArea';
                case 'Fixed Area', userInterface_.seedStrategy = 'FixedArea';
            end
            userInterface_.seedDensCtrl = app.SeedDensityCtrlEditField_2.Value;
            userInterface_.selectedPrincipalStressField = [];
            if app.MajorCheckBox.Value
                userInterface_.selectedPrincipalStressField(1,end+1) = 1;
            end
            if app.MediumCheckBox.Value
                userInterface_.selectedPrincipalStressField(1,end+1) = 2;
            end            
            if app.MinorCheckBox.Value
                userInterface_.selectedPrincipalStressField(1,end+1) = 3;
            end
            if isempty(userInterface_.selectedPrincipalStressField)
                error('There is no Principal Stress Field Selected!');
            end
            if app.StepsizeScalingFacdelta_sEditField.Value<0 || app.StepsizeScalingFacdelta_sEditField.Value>1
                warning('The Stepsize Scaling Factor must be larger than 0 and less than 1!');
                app.StepsizeScalingFacdelta_sEditField.Value = 1.0;
            end
            integratingStepScalingFac_ = app.StepsizeScalingFacdelta_sEditField.Value;
            
            userInterface_.mergingOpt = app.MergingCheckBox.Value;
            userInterface_.snappingOpt = app.SnappingCheckBox.Value;
            userInterface_.maxAngleDevi = app.PermittedMaxTangentDeviEditField.Value;
            userInterface_.multiMergingThresholds(1) = [app.MergingThresholdScalingFacMajorEditField.Value];
            userInterface_.multiMergingThresholds(2) = [app.MergingThresholdScalingFacMediumEditField.Value];
            userInterface_.multiMergingThresholds(3) = [app.MergingThresholdScalingFacMinorEditField.Value];
            switch app.IntegratingSchemeDropDown.Value
                case 'Euler'
                    userInterface_.traceAlgorithm = 'Euler';
                case 'Runge-Kutta (2)'
                    userInterface_.traceAlgorithm = 'RK2';
                case 'Runge-Kutta (4)'
                    userInterface_.traceAlgorithm = 'RK4';
            end
            RunMission(userInterface_);
            
            %%Reset Interaction Options
            app.MajorLoDsDropDown.Value = 'Geo';
            app.MajorPSLGeometryDropDown.Value = 'TUBE';
            app.MajorResCtrlSlider.Value = 0;
            app.MediumLoDsDropDown.Value = 'Geo';
            app.MediumPSLGeometryDropDown.Value = 'TUBE';
            app.MediumResCtrlSlider.Value = 0;
            app.MinorLoDsDropDown.Value = 'Geo';
            app.MinorPSLGeometryDropDown.Value = 'TUBE';
            app.MinorResCtrlSlider.Value = 0;
            app.ColoringDropDown.Value = 'None';
            
            %%Vis.
            PSLGeometryScalingFacEditFieldValueChanged(app, event);            
        end

        % Value changed function: PSLGeometryScalingFacEditField
        function PSLGeometryScalingFacEditFieldValueChanged(app, event)
           lw = app.PSLGeometryScalingFacEditField.Value;
           CreatePSLsGeometry(lw);
           PSLsMenuSelected(app, event);
        end

        % Menu selected function: ProblemDescriptionMenu
        function ProblemDescriptionMenuSelected(app, event)
            global runEnvironment_;
            global axHandle_;
            global handleSilhoutte_;
            global handleBC_;
            if isempty(runEnvironment_), return; end
            [az, el] = view(axHandle_);
            cla(axHandle_);
            colorbar(axHandle_, 'off'); 
            handleSilhoutte_ = DrawSilhouette(axHandle_);
            handleBC_ = DrawBoundaryCondition(axHandle_);
            view(axHandle_, az, el);            
            camlight(axHandle_, 'headlight','infinite');             
        end

        % Menu selected function: MeshMenu
        function MeshMenuSelected(app, event)
            global runEnvironment_;
            global axHandle_;
            if isempty(runEnvironment_), return; end
            [az, el] = view(axHandle_);
            cla(axHandle_);
            colorbar(axHandle_, 'off'); 
            DrawMesh(axHandle_);
            view(axHandle_, az, el);            
        end

        % Menu selected function: vonMisesStressMenu
        function vonMisesStressMenuSelected(app, event)
            global runEnvironment_;
            global axHandle_;
            if isempty(runEnvironment_), return; end
            [az, el] = view(axHandle_);
            cla(axHandle_);
            colorbar(axHandle_, 'off'); 
            DrawScalarStressField('Sigma_vM', axHandle_);
            view(axHandle_, az, el);            
        end

        % Menu selected function: xxMenu
        function xxMenuSelected(app, event)
            global runEnvironment_;
            global axHandle_;
            if isempty(runEnvironment_), return; end
            [az, el] = view(axHandle_);
            cla(axHandle_);
            colorbar(axHandle_, 'off'); 
            DrawScalarStressField('Sigma_xx', axHandle_);
            view(axHandle_, az, el);             
        end

        % Menu selected function: yyMenu
        function yyMenuSelected(app, event)
            global runEnvironment_;
            global axHandle_;
            if isempty(runEnvironment_), return; end
            [az, el] = view(axHandle_);
            cla(axHandle_);
            colorbar(axHandle_, 'off'); 
            DrawScalarStressField('Sigma_yy', axHandle_);
            view(axHandle_, az, el);             
        end

        % Menu selected function: zzMenu
        function zzMenuSelected(app, event)
            global runEnvironment_;
            global axHandle_;
            if isempty(runEnvironment_), return; end
            [az, el] = view(axHandle_);
            cla(axHandle_);
            colorbar(axHandle_, 'off'); 
            DrawScalarStressField('Sigma_zz', axHandle_);
            view(axHandle_, az, el);             
        end

        % Menu selected function: yzMenu
        function yzMenuSelected(app, event)
            global runEnvironment_;
            global axHandle_;
            if isempty(runEnvironment_), return; end
            [az, el] = view(axHandle_);
            cla(axHandle_);
            colorbar(axHandle_, 'off'); 
            DrawScalarStressField('Sigma_yz', axHandle_);
            view(axHandle_, az, el);             
        end

        % Menu selected function: zxMenu
        function zxMenuSelected(app, event)
            global runEnvironment_;
            global axHandle_;
            if isempty(runEnvironment_), return; end
            [az, el] = view(axHandle_);
            cla(axHandle_);
            colorbar(axHandle_, 'off'); 
            DrawScalarStressField('Sigma_zx', axHandle_);
            view(axHandle_, az, el);             
        end

        % Menu selected function: xyMenu
        function xyMenuSelected(app, event)
            global runEnvironment_;
            global axHandle_;
            if isempty(runEnvironment_), return; end
            [az, el] = view(axHandle_);
            cla(axHandle_);
            colorbar(axHandle_, 'off'); 
            DrawScalarStressField('Sigma_xy', axHandle_);
            view(axHandle_, az, el);             
        end

        % Menu selected function: MajorMenu
        function MajorMenuSelected(app, event)
            global runEnvironment_;
            global axHandle_;
            if isempty(runEnvironment_), return; end
            [az, el] = view(axHandle_);
            cla(axHandle_);
            colorbar(axHandle_, 'off'); 
            DrawScalarStressField('Sigma_1', axHandle_);
            view(axHandle_, az, el);             
        end

        % Menu selected function: MediumMenu
        function MediumMenuSelected(app, event)
            global runEnvironment_;
            global axHandle_;
            if isempty(runEnvironment_), return; end
            [az, el] = view(axHandle_);
            cla(axHandle_);
            colorbar(axHandle_, 'off'); 
            DrawScalarStressField('Sigma_2', axHandle_);
            view(axHandle_, az, el);             
        end

        % Menu selected function: MinorMenu
        function MinorMenuSelected(app, event)
            global runEnvironment_;
            global axHandle_;
            if isempty(runEnvironment_), return; end
            [az, el] = view(axHandle_);
            cla(axHandle_);
            colorbar(axHandle_, 'off'); 
            DrawScalarStressField('Sigma_3', axHandle_);
            view(axHandle_, az, el);             
        end

        % Menu selected function: SeedsMenu
        function SeedsMenuSelected(app, event)
            global runEnvironment_;
            global axHandle_;
            if isempty(runEnvironment_), return; end
            [az, el] = view(axHandle_);
            cla(axHandle_);
            colorbar(axHandle_, 'off'); 
            DrawSeedPoints(app.PSLGeometryScalingFacEditField.Value, 'inputSeeds', axHandle_);
            view(axHandle_, az, el);              
        end

        % Callback function: ColoringDropDown, MajorLoDsDropDown, 
        % MajorPSLGeometryDropDown, MajorResCtrlSlider, 
        % MediumLoDsDropDown, MediumPSLGeometryDropDown, 
        % MediumResCtrlSlider, MinorLoDsDropDown, 
        % MinorPSLGeometryDropDown, MinorResCtrlSlider, PSLsMenu, 
        % PermittedMinLengthofVisiblePSLEditField
        function PSLsMenuSelected(app, event)
            global runEnvironment_;
            global axHandle_;
	        global majorPSLpool_; 
	        global mediumPSLpool_; 
	        global minorPSLpool_;            
            if isempty(runEnvironment_), return; end
			optMajorPSLs = 0;
			if isempty(majorPSLpool_), optMajorPSLs = 1;
			else
				if isempty(majorPSLpool_(1).phyCoordList), optMajorPSLs = 1; end
			end
			optMediumPSLs = 0;
			if isempty(mediumPSLpool_), optMediumPSLs = 1;
			else
				if isempty(mediumPSLpool_(1).phyCoordList), optMediumPSLs = 1; end
			end
			optMinorPSLs = 0;
			if isempty(minorPSLpool_), optMinorPSLs = 1;
			else
				if isempty(minorPSLpool_(1).phyCoordList), optMinorPSLs = 1; end
			end			
            if optMajorPSLs && optMediumPSLs && optMinorPSLs, return; end
            [az, el] = view(axHandle_);
            cla(axHandle_);
            colorbar(axHandle_, 'off');            
            [imOpt, imVal, pslGeo, stressComponentOpt, miniPSLength] = GatherPSLvisCtrls(app);
            DrawPSLs_viaGUI(imOpt, imVal, pslGeo, stressComponentOpt, miniPSLength, axHandle_);           
            view(axHandle_, az, el);            
        end

        % Value changed function: CameraZoominoutEditField
        function CameraZoominoutEditFieldValueChanged(app, event)
            value = app.CameraZoominoutEditField.Value;
            global axHandle_;
            if value<=0
                warning('Wrong Input! Must be Positive!');
                app.CameraZoominoutEditField.Value = 1.0;
            end
            camzoom(axHandle_, app.CameraZoominoutEditField.Value);
            app.CameraZoominoutEditField.Value = 1.0;
        end

        % Menu selected function: ExportPSLsforLineVisdatMenu
        function ExportPSLsforLineVisdatMenuSelected(app, event)
            global runEnvironment_;
	        global majorPSLpool_; 
	        global mediumPSLpool_; 
	        global minorPSLpool_;
            if isempty(runEnvironment_), return; end
            if isempty(majorPSLpool_) && isempty(mediumPSLpool_) && isempty(minorPSLpool_), return; end
            [fileName, dataPath] = uiputfile('*.dat', 'Select a Path to Write the PSL Data for LineVis');
            pslDataNameOutput = strcat(dataPath, fileName);
            ExportResult(pslDataNameOutput);             
        end

        % Menu selected function: ImportStressFieldMenu
        function ImportStressFieldMenuSelected(app, event)
            addpath('../backend');
            global axHandle_;
            global userInterface_;
            global runEnvironment_;
            clc;
            runEnvironment_ = 0;
            [fileName, dataPath] = uigetfile({'*.TSVcarti'; '*.TSV'}, 'Select a Stress Field File to Open');
            %[fileName, dataPath] = uigetfile('*.*', 'Select a Stress Field File to Open');
            if isnumeric(fileName) || isnumeric(dataPath), return; end
            userInterface_ = InterfaceStruct();
            userInterface_.fileName = strcat(dataPath,fileName);
            GlobalVariables;
            disp('Loading Dataset....');
            ImportStressFields(userInterface_.fileName);
            dataName_ = userInterface_.fileName;
            axHandle_ = app.UIAxes;
            view(axHandle_, 3);
            ProblemDescriptionMenuSelected(app, event);            
            runEnvironment_ = 1;            
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {778, 778};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {331, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 1150 778];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {331, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create SubmitButton
            app.SubmitButton = uibutton(app.LeftPanel, 'push');
            app.SubmitButton.ButtonPushedFcn = createCallbackFcn(app, @SubmitButtonPushed, true);
            app.SubmitButton.BackgroundColor = [1 0.4118 0.1608];
            app.SubmitButton.Position = [124 665 100 22];
            app.SubmitButton.Text = 'Submit';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.LeftPanel);
            app.TabGroup.Position = [6 15 320 627];

            % Create SettingsTab
            app.SettingsTab = uitab(app.TabGroup);
            app.SettingsTab.Title = 'Settings';

            % Create SnappingCheckBox
            app.SnappingCheckBox = uicheckbox(app.SettingsTab);
            app.SnappingCheckBox.Text = 'Snapping';
            app.SnappingCheckBox.Position = [224 438 73 22];

            % Create MergingCheckBox
            app.MergingCheckBox = uicheckbox(app.SettingsTab);
            app.MergingCheckBox.Text = 'Merging';
            app.MergingCheckBox.Position = [27 438 65 22];
            app.MergingCheckBox.Value = true;

            % Create MinorCheckBox
            app.MinorCheckBox = uicheckbox(app.SettingsTab);
            app.MinorCheckBox.Text = 'Minor';
            app.MinorCheckBox.Position = [254 567 52 22];
            app.MinorCheckBox.Value = true;

            % Create MediumCheckBox
            app.MediumCheckBox = uicheckbox(app.SettingsTab);
            app.MediumCheckBox.Text = 'Medium';
            app.MediumCheckBox.Position = [132 567 65 22];

            % Create MajorCheckBox
            app.MajorCheckBox = uicheckbox(app.SettingsTab);
            app.MajorCheckBox.Text = 'Major';
            app.MajorCheckBox.Position = [23 567 52 22];
            app.MajorCheckBox.Value = true;

            % Create SeedingStrategyDropDown
            app.SeedingStrategyDropDown = uidropdown(app.SettingsTab);
            app.SeedingStrategyDropDown.Items = {'Volume', 'Surface', 'Loading Area', 'Fixed Area'};
            app.SeedingStrategyDropDown.Position = [157 525 134 22];
            app.SeedingStrategyDropDown.Value = 'Volume';

            % Create SeedingStrategyDropDownLabel
            app.SeedingStrategyDropDownLabel = uilabel(app.SettingsTab);
            app.SeedingStrategyDropDownLabel.HorizontalAlignment = 'right';
            app.SeedingStrategyDropDownLabel.Position = [44 525 98 22];
            app.SeedingStrategyDropDownLabel.Text = 'Seeding Strategy';

            % Create PermittedMaxTangentDeviEditFieldLabel
            app.PermittedMaxTangentDeviEditFieldLabel = uilabel(app.SettingsTab);
            app.PermittedMaxTangentDeviEditFieldLabel.HorizontalAlignment = 'right';
            app.PermittedMaxTangentDeviEditFieldLabel.Position = [87 348 163 22];
            app.PermittedMaxTangentDeviEditFieldLabel.Text = 'Permitted Max. Tangent Devi.';

            % Create PermittedMaxTangentDeviEditField
            app.PermittedMaxTangentDeviEditField = uieditfield(app.SettingsTab, 'numeric');
            app.PermittedMaxTangentDeviEditField.Position = [262 348 40 22];
            app.PermittedMaxTangentDeviEditField.Value = 20;

            % Create IntegratingSchemeDropDownLabel
            app.IntegratingSchemeDropDownLabel = uilabel(app.SettingsTab);
            app.IntegratingSchemeDropDownLabel.HorizontalAlignment = 'right';
            app.IntegratingSchemeDropDownLabel.Position = [12 398 110 22];
            app.IntegratingSchemeDropDownLabel.Text = 'Integrating Scheme';

            % Create IntegratingSchemeDropDown
            app.IntegratingSchemeDropDown = uidropdown(app.SettingsTab);
            app.IntegratingSchemeDropDown.Items = {'Euler', 'Runge-Kutta (2)', 'Runge-Kutta (4)'};
            app.IntegratingSchemeDropDown.Position = [137 398 158 22];
            app.IntegratingSchemeDropDown.Value = 'Euler';

            % Create MergingThresholdScalingFacMajorEditFieldLabel
            app.MergingThresholdScalingFacMajorEditFieldLabel = uilabel(app.SettingsTab);
            app.MergingThresholdScalingFacMajorEditFieldLabel.BackgroundColor = [0.651 0.651 0.651];
            app.MergingThresholdScalingFacMajorEditFieldLabel.HorizontalAlignment = 'right';
            app.MergingThresholdScalingFacMajorEditFieldLabel.Position = [26 308 225 22];
            app.MergingThresholdScalingFacMajorEditFieldLabel.Text = 'Merging Threshold Scaling. Fac. (Major) ';

            % Create MergingThresholdScalingFacMajorEditField
            app.MergingThresholdScalingFacMajorEditField = uieditfield(app.SettingsTab, 'numeric');
            app.MergingThresholdScalingFacMajorEditField.BackgroundColor = [0.651 0.651 0.651];
            app.MergingThresholdScalingFacMajorEditField.Position = [263 308 40 22];
            app.MergingThresholdScalingFacMajorEditField.Value = 1;

            % Create MergingThresholdScalingFacMediumEditFieldLabel
            app.MergingThresholdScalingFacMediumEditFieldLabel = uilabel(app.SettingsTab);
            app.MergingThresholdScalingFacMediumEditFieldLabel.BackgroundColor = [0.651 0.651 0.651];
            app.MergingThresholdScalingFacMediumEditFieldLabel.HorizontalAlignment = 'right';
            app.MergingThresholdScalingFacMediumEditFieldLabel.Position = [12 266 238 22];
            app.MergingThresholdScalingFacMediumEditFieldLabel.Text = 'Merging Threshold Scaling. Fac. (Medium) ';

            % Create MergingThresholdScalingFacMediumEditField
            app.MergingThresholdScalingFacMediumEditField = uieditfield(app.SettingsTab, 'numeric');
            app.MergingThresholdScalingFacMediumEditField.BackgroundColor = [0.651 0.651 0.651];
            app.MergingThresholdScalingFacMediumEditField.Position = [262 266 40 22];
            app.MergingThresholdScalingFacMediumEditField.Value = 1;

            % Create MergingThresholdScalingFacMinorEditFieldLabel
            app.MergingThresholdScalingFacMinorEditFieldLabel = uilabel(app.SettingsTab);
            app.MergingThresholdScalingFacMinorEditFieldLabel.BackgroundColor = [0.651 0.651 0.651];
            app.MergingThresholdScalingFacMinorEditFieldLabel.HorizontalAlignment = 'right';
            app.MergingThresholdScalingFacMinorEditFieldLabel.Position = [26 225 225 22];
            app.MergingThresholdScalingFacMinorEditFieldLabel.Text = 'Merging Threshold Scaling. Fac. (Minor) ';

            % Create MergingThresholdScalingFacMinorEditField
            app.MergingThresholdScalingFacMinorEditField = uieditfield(app.SettingsTab, 'numeric');
            app.MergingThresholdScalingFacMinorEditField.BackgroundColor = [0.651 0.651 0.651];
            app.MergingThresholdScalingFacMinorEditField.Position = [263 225 40 22];
            app.MergingThresholdScalingFacMinorEditField.Value = 1;

            % Create PSLGeometryScalingFacLabel
            app.PSLGeometryScalingFacLabel = uilabel(app.SettingsTab);
            app.PSLGeometryScalingFacLabel.HorizontalAlignment = 'right';
            app.PSLGeometryScalingFacLabel.Position = [98 184 154 22];
            app.PSLGeometryScalingFacLabel.Text = 'PSL Geometry Scaling Fac.';

            % Create PSLGeometryScalingFacEditField
            app.PSLGeometryScalingFacEditField = uieditfield(app.SettingsTab, 'numeric');
            app.PSLGeometryScalingFacEditField.ValueChangedFcn = createCallbackFcn(app, @PSLGeometryScalingFacEditFieldValueChanged, true);
            app.PSLGeometryScalingFacEditField.Position = [264 184 40 22];
            app.PSLGeometryScalingFacEditField.Value = 0.5;

            % Create PermittedMinLengthofVisiblePSLEditFieldLabel
            app.PermittedMinLengthofVisiblePSLEditFieldLabel = uilabel(app.SettingsTab);
            app.PermittedMinLengthofVisiblePSLEditFieldLabel.HorizontalAlignment = 'right';
            app.PermittedMinLengthofVisiblePSLEditFieldLabel.Position = [50 143 202 22];
            app.PermittedMinLengthofVisiblePSLEditFieldLabel.Text = 'Permitted Min. Length of Visible PSL';

            % Create PermittedMinLengthofVisiblePSLEditField
            app.PermittedMinLengthofVisiblePSLEditField = uieditfield(app.SettingsTab, 'numeric');
            app.PermittedMinLengthofVisiblePSLEditField.ValueChangedFcn = createCallbackFcn(app, @PSLsMenuSelected, true);
            app.PermittedMinLengthofVisiblePSLEditField.Position = [264 143 40 22];
            app.PermittedMinLengthofVisiblePSLEditField.Value = 10;

            % Create VoxelizedPSLThicknessExportEditFieldLabel
            app.VoxelizedPSLThicknessExportEditFieldLabel = uilabel(app.SettingsTab);
            app.VoxelizedPSLThicknessExportEditFieldLabel.BackgroundColor = [0.502 0.502 0.502];
            app.VoxelizedPSLThicknessExportEditFieldLabel.HorizontalAlignment = 'right';
            app.VoxelizedPSLThicknessExportEditFieldLabel.Position = [65 51 186 22];
            app.VoxelizedPSLThicknessExportEditFieldLabel.Text = 'Voxelized PSL Thickness (Export)';

            % Create VoxelizedPSLThicknessExportEditField
            app.VoxelizedPSLThicknessExportEditField = uieditfield(app.SettingsTab, 'numeric');
            app.VoxelizedPSLThicknessExportEditField.BackgroundColor = [0.502 0.502 0.502];
            app.VoxelizedPSLThicknessExportEditField.Position = [263 51 40 22];
            app.VoxelizedPSLThicknessExportEditField.Value = 1;

            % Create BoundaryThicknessExportEditFieldLabel
            app.BoundaryThicknessExportEditFieldLabel = uilabel(app.SettingsTab);
            app.BoundaryThicknessExportEditFieldLabel.BackgroundColor = [0.502 0.502 0.502];
            app.BoundaryThicknessExportEditFieldLabel.HorizontalAlignment = 'right';
            app.BoundaryThicknessExportEditFieldLabel.Position = [90 15 160 22];
            app.BoundaryThicknessExportEditFieldLabel.Text = 'Boundary Thickness (Export)';

            % Create BoundaryThicknessExportEditField
            app.BoundaryThicknessExportEditField = uieditfield(app.SettingsTab, 'numeric');
            app.BoundaryThicknessExportEditField.BackgroundColor = [0.502 0.502 0.502];
            app.BoundaryThicknessExportEditField.Position = [262 15 40 22];
            app.BoundaryThicknessExportEditField.Value = 1;

            % Create SeedDensityCtrlEditField_2Label
            app.SeedDensityCtrlEditField_2Label = uilabel(app.SettingsTab);
            app.SeedDensityCtrlEditField_2Label.HorizontalAlignment = 'right';
            app.SeedDensityCtrlEditField_2Label.Position = [44 479 102 22];
            app.SeedDensityCtrlEditField_2Label.Text = 'Seed Density Ctrl.';

            % Create SeedDensityCtrlEditField_2
            app.SeedDensityCtrlEditField_2 = uieditfield(app.SettingsTab, 'numeric');
            app.SeedDensityCtrlEditField_2.Position = [161 479 100 22];
            app.SeedDensityCtrlEditField_2.Value = 5;

            % Create StepsizeScalingFacdelta_sEditFieldLabel
            app.StepsizeScalingFacdelta_sEditFieldLabel = uilabel(app.SettingsTab);
            app.StepsizeScalingFacdelta_sEditFieldLabel.HorizontalAlignment = 'right';
            app.StepsizeScalingFacdelta_sEditFieldLabel.Position = [80 102 171 22];
            app.StepsizeScalingFacdelta_sEditFieldLabel.Text = 'Stepsize Scaling Fac. (delta_s)';

            % Create StepsizeScalingFacdelta_sEditField
            app.StepsizeScalingFacdelta_sEditField = uieditfield(app.SettingsTab, 'numeric');
            app.StepsizeScalingFacdelta_sEditField.Position = [263 102 40 22];
            app.StepsizeScalingFacdelta_sEditField.Value = 1;

            % Create InteractionsTab
            app.InteractionsTab = uitab(app.TabGroup);
            app.InteractionsTab.Title = 'Interactions';

            % Create MajorLoDsDropDownLabel
            app.MajorLoDsDropDownLabel = uilabel(app.InteractionsTab);
            app.MajorLoDsDropDownLabel.HorizontalAlignment = 'right';
            app.MajorLoDsDropDownLabel.Position = [92 567 67 22];
            app.MajorLoDsDropDownLabel.Text = 'Major LoDs';

            % Create MajorLoDsDropDown
            app.MajorLoDsDropDown = uidropdown(app.InteractionsTab);
            app.MajorLoDsDropDown.Items = {'Geo', 'PS', 'vM', 'Length'};
            app.MajorLoDsDropDown.ValueChangedFcn = createCallbackFcn(app, @PSLsMenuSelected, true);
            app.MajorLoDsDropDown.Position = [174 567 115 22];
            app.MajorLoDsDropDown.Value = 'Geo';

            % Create MajorPSLGeometryDropDownLabel
            app.MajorPSLGeometryDropDownLabel = uilabel(app.InteractionsTab);
            app.MajorPSLGeometryDropDownLabel.HorizontalAlignment = 'right';
            app.MajorPSLGeometryDropDownLabel.Position = [91 533 118 22];
            app.MajorPSLGeometryDropDownLabel.Text = 'Major PSL Geometry';

            % Create MajorPSLGeometryDropDown
            app.MajorPSLGeometryDropDown = uidropdown(app.InteractionsTab);
            app.MajorPSLGeometryDropDown.Items = {'TUBE', 'RIBBON'};
            app.MajorPSLGeometryDropDown.ValueChangedFcn = createCallbackFcn(app, @PSLsMenuSelected, true);
            app.MajorPSLGeometryDropDown.Position = [224 533 65 22];
            app.MajorPSLGeometryDropDown.Value = 'TUBE';

            % Create MajorResCtrlSliderLabel
            app.MajorResCtrlSliderLabel = uilabel(app.InteractionsTab);
            app.MajorResCtrlSliderLabel.HorizontalAlignment = 'right';
            app.MajorResCtrlSliderLabel.Position = [18 500 86 22];
            app.MajorResCtrlSliderLabel.Text = 'Major Res. Ctrl';

            % Create MajorResCtrlSlider
            app.MajorResCtrlSlider = uislider(app.InteractionsTab);
            app.MajorResCtrlSlider.Limits = [0 1.2];
            app.MajorResCtrlSlider.ValueChangedFcn = createCallbackFcn(app, @PSLsMenuSelected, true);
            app.MajorResCtrlSlider.Position = [125 509 158 3];

            % Create MediumResCtrlSliderLabel
            app.MediumResCtrlSliderLabel = uilabel(app.InteractionsTab);
            app.MediumResCtrlSliderLabel.HorizontalAlignment = 'right';
            app.MediumResCtrlSliderLabel.Position = [7 359 98 22];
            app.MediumResCtrlSliderLabel.Text = 'Medium Res. Ctrl';

            % Create MediumResCtrlSlider
            app.MediumResCtrlSlider = uislider(app.InteractionsTab);
            app.MediumResCtrlSlider.Limits = [0 1.2];
            app.MediumResCtrlSlider.ValueChangedFcn = createCallbackFcn(app, @PSLsMenuSelected, true);
            app.MediumResCtrlSlider.Position = [126 368 157 3];

            % Create MediumLoDsDropDownLabel
            app.MediumLoDsDropDownLabel = uilabel(app.InteractionsTab);
            app.MediumLoDsDropDownLabel.HorizontalAlignment = 'right';
            app.MediumLoDsDropDownLabel.Position = [79 428 80 22];
            app.MediumLoDsDropDownLabel.Text = 'Medium LoDs';

            % Create MediumLoDsDropDown
            app.MediumLoDsDropDown = uidropdown(app.InteractionsTab);
            app.MediumLoDsDropDown.Items = {'Geo', 'PS', 'vM', 'Length'};
            app.MediumLoDsDropDown.ValueChangedFcn = createCallbackFcn(app, @PSLsMenuSelected, true);
            app.MediumLoDsDropDown.Position = [174 428 115 22];
            app.MediumLoDsDropDown.Value = 'Geo';

            % Create MediumPSLGeometryDropDownLabel
            app.MediumPSLGeometryDropDownLabel = uilabel(app.InteractionsTab);
            app.MediumPSLGeometryDropDownLabel.HorizontalAlignment = 'right';
            app.MediumPSLGeometryDropDownLabel.Position = [78 389 130 22];
            app.MediumPSLGeometryDropDownLabel.Text = 'Medium PSL Geometry';

            % Create MediumPSLGeometryDropDown
            app.MediumPSLGeometryDropDown = uidropdown(app.InteractionsTab);
            app.MediumPSLGeometryDropDown.Items = {'TUBE', 'RIBBON'};
            app.MediumPSLGeometryDropDown.ValueChangedFcn = createCallbackFcn(app, @PSLsMenuSelected, true);
            app.MediumPSLGeometryDropDown.Position = [223 389 65 22];
            app.MediumPSLGeometryDropDown.Value = 'TUBE';

            % Create MinorResCtrlSliderLabel
            app.MinorResCtrlSliderLabel = uilabel(app.InteractionsTab);
            app.MinorResCtrlSliderLabel.HorizontalAlignment = 'right';
            app.MinorResCtrlSliderLabel.Position = [17 214 86 22];
            app.MinorResCtrlSliderLabel.Text = 'Minor Res. Ctrl';

            % Create MinorResCtrlSlider
            app.MinorResCtrlSlider = uislider(app.InteractionsTab);
            app.MinorResCtrlSlider.Limits = [0 1.2];
            app.MinorResCtrlSlider.ValueChangedFcn = createCallbackFcn(app, @PSLsMenuSelected, true);
            app.MinorResCtrlSlider.Position = [124 223 158 3];

            % Create MinorLoDsDropDownLabel
            app.MinorLoDsDropDownLabel = uilabel(app.InteractionsTab);
            app.MinorLoDsDropDownLabel.HorizontalAlignment = 'right';
            app.MinorLoDsDropDownLabel.Position = [87 287 67 22];
            app.MinorLoDsDropDownLabel.Text = 'Minor LoDs';

            % Create MinorLoDsDropDown
            app.MinorLoDsDropDown = uidropdown(app.InteractionsTab);
            app.MinorLoDsDropDown.Items = {'Geo', 'PS', 'vM', 'Length'};
            app.MinorLoDsDropDown.ValueChangedFcn = createCallbackFcn(app, @PSLsMenuSelected, true);
            app.MinorLoDsDropDown.Position = [169 287 115 22];
            app.MinorLoDsDropDown.Value = 'Geo';

            % Create MinorPSLGeometryDropDownLabel
            app.MinorPSLGeometryDropDownLabel = uilabel(app.InteractionsTab);
            app.MinorPSLGeometryDropDownLabel.HorizontalAlignment = 'right';
            app.MinorPSLGeometryDropDownLabel.Position = [85 246 118 22];
            app.MinorPSLGeometryDropDownLabel.Text = 'Minor PSL Geometry';

            % Create MinorPSLGeometryDropDown
            app.MinorPSLGeometryDropDown = uidropdown(app.InteractionsTab);
            app.MinorPSLGeometryDropDown.Items = {'TUBE', 'RIBBON'};
            app.MinorPSLGeometryDropDown.ValueChangedFcn = createCallbackFcn(app, @PSLsMenuSelected, true);
            app.MinorPSLGeometryDropDown.Position = [218 246 65 22];
            app.MinorPSLGeometryDropDown.Value = 'TUBE';

            % Create ColoringDropDownLabel
            app.ColoringDropDownLabel = uilabel(app.InteractionsTab);
            app.ColoringDropDownLabel.HorizontalAlignment = 'right';
            app.ColoringDropDownLabel.Position = [74 143 50 22];
            app.ColoringDropDownLabel.Text = 'Coloring';

            % Create ColoringDropDown
            app.ColoringDropDown = uidropdown(app.InteractionsTab);
            app.ColoringDropDown.Items = {'None', 'Principal Stress', 'Mises Stress', 'Normal Stress (xx)', 'Normal Stress (yy)', 'Normal Stress (zz)', 'Shear Stress (yz)', 'Shear Stress (zx)', 'Shear Stress (xy)'};
            app.ColoringDropDown.ValueChangedFcn = createCallbackFcn(app, @PSLsMenuSelected, true);
            app.ColoringDropDown.Position = [139 143 152 22];
            app.ColoringDropDown.Value = 'None';

            % Create CameraZoominoutEditFieldLabel
            app.CameraZoominoutEditFieldLabel = uilabel(app.InteractionsTab);
            app.CameraZoominoutEditFieldLabel.BackgroundColor = [0.651 0.651 0.651];
            app.CameraZoominoutEditFieldLabel.HorizontalAlignment = 'right';
            app.CameraZoominoutEditFieldLabel.Position = [122 93 115 22];
            app.CameraZoominoutEditFieldLabel.Text = 'Camera Zoom in/out';

            % Create CameraZoominoutEditField
            app.CameraZoominoutEditField = uieditfield(app.InteractionsTab, 'numeric');
            app.CameraZoominoutEditField.ValueChangedFcn = createCallbackFcn(app, @CameraZoominoutEditFieldValueChanged, true);
            app.CameraZoominoutEditField.BackgroundColor = [0.651 0.651 0.651];
            app.CameraZoominoutEditField.Position = [249 93 40 22];
            app.CameraZoominoutEditField.Value = 1;

            % Create LineDensityCtrlEditField_2Label
            app.LineDensityCtrlEditField_2Label = uilabel(app.LeftPanel);
            app.LineDensityCtrlEditField_2Label.HorizontalAlignment = 'right';
            app.LineDensityCtrlEditField_2Label.Position = [13 750 97 22];
            app.LineDensityCtrlEditField_2Label.Text = 'Line Density Ctrl.';

            % Create LineDensityCtrlEditField_2
            app.LineDensityCtrlEditField_2 = uieditfield(app.LeftPanel, 'numeric');
            app.LineDensityCtrlEditField_2.Position = [125 750 100 22];
            app.LineDensityCtrlEditField_2.Value = 5;

            % Create NumLevelsEditField_2Label
            app.NumLevelsEditField_2Label = uilabel(app.LeftPanel);
            app.NumLevelsEditField_2Label.HorizontalAlignment = 'right';
            app.NumLevelsEditField_2Label.Position = [37 706 72 22];
            app.NumLevelsEditField_2Label.Text = 'Num. Levels';

            % Create NumLevelsEditField_2
            app.NumLevelsEditField_2 = uieditfield(app.LeftPanel, 'numeric');
            app.NumLevelsEditField_2.Position = [124 706 100 22];
            app.NumLevelsEditField_2.Value = 1;

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create UIAxes
            app.UIAxes = uiaxes(app.RightPanel);
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.View = [-37.5 30];
            app.UIAxes.Projection = 'perspective';
            app.UIAxes.FontName = 'Times New Roman';
            app.UIAxes.FontSize = 20;
            app.UIAxes.Position = [7 66 785 685];

            % Create FilesMenu
            app.FilesMenu = uimenu(app.UIFigure);
            app.FilesMenu.Text = 'Files';

            % Create ImportStressFieldMenu
            app.ImportStressFieldMenu = uimenu(app.FilesMenu);
            app.ImportStressFieldMenu.MenuSelectedFcn = createCallbackFcn(app, @ImportStressFieldMenuSelected, true);
            app.ImportStressFieldMenu.Text = 'Import Stress Field';

            % Create ExportPSLsforLineVisdatMenu
            app.ExportPSLsforLineVisdatMenu = uimenu(app.FilesMenu);
            app.ExportPSLsforLineVisdatMenu.MenuSelectedFcn = createCallbackFcn(app, @ExportPSLsforLineVisdatMenuSelected, true);
            app.ExportPSLsforLineVisdatMenu.Text = 'Export PSLs for "LineVis" (*.dat)';

            % Create MediaMenu
            app.MediaMenu = uimenu(app.UIFigure);
            app.MediaMenu.Text = 'Media';

            % Create ProblemDescriptionMenu
            app.ProblemDescriptionMenu = uimenu(app.MediaMenu);
            app.ProblemDescriptionMenu.MenuSelectedFcn = createCallbackFcn(app, @ProblemDescriptionMenuSelected, true);
            app.ProblemDescriptionMenu.Text = 'Problem Description';

            % Create MeshMenu
            app.MeshMenu = uimenu(app.MediaMenu);
            app.MeshMenu.MenuSelectedFcn = createCallbackFcn(app, @MeshMenuSelected, true);
            app.MeshMenu.Text = 'Mesh';

            % Create ScalarStressFieldMenu
            app.ScalarStressFieldMenu = uimenu(app.MediaMenu);
            app.ScalarStressFieldMenu.Text = 'Scalar Stress Field';

            % Create vonMisesStressMenu
            app.vonMisesStressMenu = uimenu(app.ScalarStressFieldMenu);
            app.vonMisesStressMenu.MenuSelectedFcn = createCallbackFcn(app, @vonMisesStressMenuSelected, true);
            app.vonMisesStressMenu.Text = 'von Mises Stress';

            % Create NormalStressMenu
            app.NormalStressMenu = uimenu(app.ScalarStressFieldMenu);
            app.NormalStressMenu.Text = 'Normal Stress';

            % Create xxMenu
            app.xxMenu = uimenu(app.NormalStressMenu);
            app.xxMenu.MenuSelectedFcn = createCallbackFcn(app, @xxMenuSelected, true);
            app.xxMenu.Text = 'xx';

            % Create yyMenu
            app.yyMenu = uimenu(app.NormalStressMenu);
            app.yyMenu.MenuSelectedFcn = createCallbackFcn(app, @yyMenuSelected, true);
            app.yyMenu.Text = 'yy';

            % Create zzMenu
            app.zzMenu = uimenu(app.NormalStressMenu);
            app.zzMenu.MenuSelectedFcn = createCallbackFcn(app, @zzMenuSelected, true);
            app.zzMenu.Text = 'zz';

            % Create ShearStressMenu
            app.ShearStressMenu = uimenu(app.ScalarStressFieldMenu);
            app.ShearStressMenu.Text = 'Shear Stress';

            % Create yzMenu
            app.yzMenu = uimenu(app.ShearStressMenu);
            app.yzMenu.MenuSelectedFcn = createCallbackFcn(app, @yzMenuSelected, true);
            app.yzMenu.Text = 'yz';

            % Create zxMenu
            app.zxMenu = uimenu(app.ShearStressMenu);
            app.zxMenu.MenuSelectedFcn = createCallbackFcn(app, @zxMenuSelected, true);
            app.zxMenu.Text = 'zx';

            % Create xyMenu
            app.xyMenu = uimenu(app.ShearStressMenu);
            app.xyMenu.MenuSelectedFcn = createCallbackFcn(app, @xyMenuSelected, true);
            app.xyMenu.Text = 'xy';

            % Create PrincipalStressMenu
            app.PrincipalStressMenu = uimenu(app.ScalarStressFieldMenu);
            app.PrincipalStressMenu.Text = 'Principal Stress';

            % Create MajorMenu
            app.MajorMenu = uimenu(app.PrincipalStressMenu);
            app.MajorMenu.MenuSelectedFcn = createCallbackFcn(app, @MajorMenuSelected, true);
            app.MajorMenu.Text = 'Major';

            % Create MediumMenu
            app.MediumMenu = uimenu(app.PrincipalStressMenu);
            app.MediumMenu.MenuSelectedFcn = createCallbackFcn(app, @MediumMenuSelected, true);
            app.MediumMenu.Text = 'Medium';

            % Create MinorMenu
            app.MinorMenu = uimenu(app.PrincipalStressMenu);
            app.MinorMenu.MenuSelectedFcn = createCallbackFcn(app, @MinorMenuSelected, true);
            app.MinorMenu.Text = 'Minor';

            % Create SeedsMenu
            app.SeedsMenu = uimenu(app.MediaMenu);
            app.SeedsMenu.MenuSelectedFcn = createCallbackFcn(app, @SeedsMenuSelected, true);
            app.SeedsMenu.Text = 'Seeds';

            % Create PSLsMenu
            app.PSLsMenu = uimenu(app.MediaMenu);
            app.PSLsMenu.MenuSelectedFcn = createCallbackFcn(app, @PSLsMenuSelected, true);
            app.PSLsMenu.Text = 'PSLs';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = TSV3D_GUI

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end