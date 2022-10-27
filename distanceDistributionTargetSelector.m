function distanceDistributionTarget = distanceDistributionTargetSelector(app)
    distanceDistributionTargetProperties_temp = drawpoint(app.AxesForPrTargeting, 'Visible', 'off');
    
    distanceDistributionTarget = distanceDistributionTargetProperties_temp.Position(1)-0.03;
    app.TargetDistanceEditField.Value = distanceDistributionTargetProperties_temp.Position(1)-0.03;
end