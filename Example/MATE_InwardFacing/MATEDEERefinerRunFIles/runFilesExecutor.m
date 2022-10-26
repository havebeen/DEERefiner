function DEERRefineJob = runFilesExecutor(app)
    
    executeTimes = length(app.runFileName);
    for executeLoopNumbers = 1:executeTimes
        fileName(executeLoopNumbers) = app.runFileName(executeLoopNumbers);
        DEERRefineJob(executeLoopNumbers) = batch(fileName(executeLoopNumbers));
        pause(3)
    end
end