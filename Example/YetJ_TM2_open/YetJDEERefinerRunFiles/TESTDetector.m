function TESTDetector
    TEST = load('variable_test.dat');
    if TEST == 1
        a = 1
%         app.EditField.Value = num2str(app.numberOfCPUCores-1);
%         app.EditField.FontColor = [0.5 0.237 0.88];
%         drawnow
    elseif TEST == 2
        a = 2
%         app.EditField.Value = num2str(app.numberOfCPUCores);
%         app.EditField.FontColor = [0 0.8 0];
%         drawnow
    end
end