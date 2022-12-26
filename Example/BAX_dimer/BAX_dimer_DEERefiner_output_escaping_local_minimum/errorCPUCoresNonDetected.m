function errorCPUCoresNonDetected
    CPUCoresNonDetectedErrorMessage = errordlg("CPU cores are not detected");
    set(CPUCoresNonDetectedErrorMessage, 'position', [100 440 450 100]);
    CPUCoresNonDetectedErrorMessageText = findobj(CPUCoresNonDetectedErrorMessage, 'Type', 'Text');
    CPUCoresNonDetectedErrorMessageText.FontSize = 25;
end