function parameterApplier(app)
    app.loopRegionString = app.FlexibleregionsEditField.Value;
    app.FlexibleregionsEditField.FontColor = [0.47 0.67 0.19];

    app.phiPsiAngle = app.phipsiangleEditField.Value;
    app.phipsiangleEditField.FontColor = [0.47 0.67 0.19];
    
    app.monteCarloIteration = app.MCiterationEditField.Value;
    app.MCiterationEditField.FontColor = [0.47 0.67 0.19];
    
    app.monteCarloTemperature = app.MCtempEditField.Value;
    app.MCtempEditField.FontColor = [0.47 0.67 0.19];
    
    app.RMSE = app.RMSEEditField.Value;
    app.RMSEEditField.FontColor = [0.47 0.67 0.19];
    
    app.maximalClashes = app.MaximalclashesEditField.Value;
    app.MaximalclashesEditField.FontColor = [0.47 0.67 0.19];
    
    app.numberOfStructure = app.NumberofstructureEditField.Value;
    app.NumberofstructureEditField.FontColor = [0.47 0.67 0.19];
    
    app.numberOfCPUCores = app.UsageofcoresEditField.Value;
    app.UsageofcoresEditField.FontColor = [0.47 0.67 0.19];
    
end