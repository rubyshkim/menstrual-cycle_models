function [E2,P4,IH,IHdel] = wright_OvarianHormones(x,lag_sol,pars)

    % get current state variable values
		RP_LH  = x(1); %current RP_LH level
		LH     = x(2); %current LH level
		RP_FSH = x(3);
		FSH    = x(4);
		RcF    = x(5); %current MsF level
		GrF    = x(6); %current GrF level
		DomF   = x(7); %current DomF level
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
        p_0       = pars(35);
        p_1       = pars(36);
        p_2       = pars(37);
        Km_Papp   = pars(38);
        mu        = pars(39);
        h_0       = pars(40);
        h_1       = pars(41);
        h_2       = pars(42);
        h_3       = pars(43);

    % Drug parameters
		p_dose      = pars(44);
		e_dose       = pars(45);
        halfLife_P   = pars(46);
		halfLife_E2  = pars(47);

    % Calculate current aux equations
		E2=e_0+e_1*GrF+e_2*DomF+e_3*Lut4+e_dose;
		P4=p_0+p_1*Lut3+p_2*Lut4+p_dose;
		IH=h_0+h_1*DomF+h_2*Lut2+h_3*Lut3;

    % Calculate delayed aux equation
      ylag1 = lag_sol(:,1);
      IHdel=h_0+h_1*ylag1(7)+h_2*ylag1(11)+h_3*ylag1(12);

    % Effective P4
      P4_eff = P4*(.5 + .5*Hill_plus(E2,Km_Papp,mu));
      P4 = P4_eff;

end
