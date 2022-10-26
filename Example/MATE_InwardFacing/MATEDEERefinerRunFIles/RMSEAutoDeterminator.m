function RMSEAutoDeterminator(app)
    distanceDistributionNumber = length(app.DEERRefineRestraintsTable.Data(:, 1));
    app.RMSE = distanceDistributionNumber*0.05+1.1;
    app.RMSEEditField.Value = app.RMSE;
    app.RMSEEditField.FontColor = [0 0 0];
end