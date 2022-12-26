function currentStatusString = currentStatusConverter(currentStatusFlag)
    if currentStatusFlag == 0
        currentStatusString = "Initializing simulations";
    elseif currentStatusFlag == 1
        currentStatusString = "Generating backbones";
    elseif currentStatusFlag == 2
        currentStatusString = "Accepting Candidate";
    elseif currentStatusFlag == 3
        currentStatusString = "RMSE criterion passed";
    elseif currentStatusFlag == 4
        currentStatusString = "Accepting Candidate";
    elseif currentStatusFlag == 5
        currentStatusString = "Clashes criterion passed";
    elseif currentStatusFlag == 6
        currentStatusString = "Generating backbones";
    elseif currentStatusFlag == 7
        currentStatusString = "Accepting Candidate";
    elseif currentStatusFlag == 8
        currentStatusString = "RMSE criterion passed";
    elseif currentStatusFlag == 9
        currentStatusString = "Accepting Candidate";
    elseif currentStatusFlag == 10
        currentStatusString = "Clashes criterion passed";
    elseif currentStatusFlag == 11
        currentStatusString = "Generating backbones";
    elseif currentStatusFlag == 12
        currentStatusString = "Accepting Candidate";
    elseif currentStatusFlag == 13
        currentStatusString = "RMSE criterion passed";
    elseif currentStatusFlag == 14
        currentStatusString = "Accepting Candidate";
    elseif currentStatusFlag == 15
        currentStatusString = "Clashes criterion passed";
    elseif currentStatusFlag == 16
        currentStatusString = "Generating backbones";
    elseif currentStatusFlag == 17
        currentStatusString = "Accepting Candidate";
    elseif currentStatusFlag == 18
        currentStatusString = "RMSE criterion passed";
    elseif currentStatusFlag == 19
        currentStatusString = "Accepting Candidate";
    elseif currentStatusFlag == 20
        currentStatusString = "Clashes criterion passed";
    elseif currentStatusFlag == 21
        currentStatusString = "Generating backbones";
    elseif currentStatusFlag == 22
        currentStatusString = "Accepting Candidate";
    elseif currentStatusFlag == 23
        currentStatusString = "Generating final answer";
    end
end