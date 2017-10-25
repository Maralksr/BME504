    function Fse = Fse(LT)
        %---------------------------
        % series elastic element (tendon)
        % input: tendon length
        % output: tendon force (0-1)
        %---------------------------
        cT_se = 27.8; 
        kT_se = 0.0047;
        LrT_se = 0.964;
        
        Fse = cT_se * kT_se * log(exp((LT - LrT_se)/kT_se)+1);
        
    end