function checkingWaitingForCPUCoresDetectedCheckBox(app)
    app.WaitingForCPUCoresDetectedCheckBox.Value = 1;  % check the box
    app.WaitingForCPUCoresDetectedCheckBox.Text = strcat(num2str(app.numberOfCPUCores), " cores detected");    % change texts in the box
    app.WaitingForCPUCoresDetectedCheckBox.FontColor = [0.47 0.67 0.19];   % make the texts green
end