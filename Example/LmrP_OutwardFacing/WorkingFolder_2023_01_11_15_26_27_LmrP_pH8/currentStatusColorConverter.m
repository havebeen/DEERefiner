function currentStatusColor = currentStatusColorConverter(currentStatusFlag)
    currentStatusColor = [(23-currentStatusFlag)/23 0 (currentStatusFlag)/23];
end