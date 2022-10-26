function StringIfCPUCoresDetected = checkIfCPUCoresDetected(numberOfCPUCores)
    BooleanIfCPUCoresDetected = isa(numberOfCPUCores,'double');
    if BooleanIfCPUCoresDetected==1
        StringIfCPUCoresDetected = "Detected";
    else
        StringIfCPUCoresDetected = "Non-detected";
    end
end