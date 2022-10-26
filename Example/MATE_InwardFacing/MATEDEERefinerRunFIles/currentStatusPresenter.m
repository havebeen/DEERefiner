function currentStatusPresenter(app)
    currentStatusFlag = currentStatusDetector(app);
    currentStatusString = currentStatusConverter(currentStatusFlag);
    currentStatusColor = currentStatusColorConverter(currentStatusFlag);
    app.EditField.Value = currentStatusString;
    app.EditField.BackgroundColor = [1 1 1];
    pause(0.01)
    app.EditField.BackgroundColor = currentStatusColor;
    passedStructureNumbers = passedStructureNumbersDetector(app);
    app.Stage1EditField.Value = passedStructureNumbers(1);
    app.Stage2EditField.Value = passedStructureNumbers(2);
    app.Stage3EditField.Value = passedStructureNumbers(3);
    app.Stage4EditField.Value = passedStructureNumbers(4);
    app.Stage5EditField.Value = passedStructureNumbers(5);
%     drawnow
end