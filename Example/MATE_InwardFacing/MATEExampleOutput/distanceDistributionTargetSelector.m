function distanceDistributionTarget = distanceDistributionTargetSelector(app)
    distanceDistributionTargetProperties_temp = drawpoint(app.AxesForPrTargeting, 'Visible', 'off');
    distanceDistributionTarget = distanceDistributionTargetProperties_temp.Position(1);
    app.TargetDistanceEditField.Value = distanceDistributionTargetProperties_temp.Position(1);
end