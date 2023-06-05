function [E2,P4,P4del,IH,IHdel] = clark_OvarianHormones(x,lag_sol,pars)

    % get current state variable values
		RP_LH  = x(1); %current RP_LH level
		LH     = x(2); %current LH level
		RP_FSH = x(3);
		FSH    = x(4);
		RcF    = x(5); %current MsF level
		SeF    = x(6); %current SeF level
		PrF    = x(7); %current PrF level
		Sc1    = x(8); %current Sc1 level
		Sc2    = x(9); %current Sc2 level
		Lut1   = x(10); %current Lut1 level
		Lut2   = x(11); %current Lut2 level
		Lut3   = x(12); %current Lut3 level
		Lut4   = x(13); %current Lut4 level

    % Get aux pars
        e_0       = pars(31);
        e_1       = pars(32);
        e_2       = pars(33);
        e_3       = pars(34);
        p_1       = pars(35);
        p_2       = pars(36);
        h_0       = pars(37);
        h_1       = pars(38);
        h_2       = pars(39);
        h_3       = pars(40);
        e_dose    = pars(41);
        p_dose    = pars(42);

    % Calculate current aux equations
		E2=e_0+e_1*SeF+e_2*PrF+e_3*Lut4 + e_dose;
		P4=p_1*Lut3+p_2*Lut4 + p_dose;
		IH=h_0+h_1*PrF+h_2*Lut3+h_3*Lut4;

    % Calculate delayed aux equations
      ylag1 = lag_sol(:,1);
      IHdel=h_0+h_1*ylag1(7)+h_2*ylag1(12)+h_3*ylag1(13);

      ylag2 = lag_sol(:,2);
      P4del= p_1*ylag2(12)+p_2*ylag2(13)+p_dose;

end
