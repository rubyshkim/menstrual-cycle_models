function [E2,P4,Inh,IHdel] = gavina_OvarianHormones(x,lag_sol,pars)
%LdeP Updated for Gavina 2022 December 29, 2022

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

    % Get aux pars - Gavina
        e0      = pars(32);
        e1      = pars(33);
        e2      = pars(34);
        e3      = pars(35);
        p1      = pars(36);
        p2      = pars(37);
        h1      = pars(38);
        h2      = pars(39);
        h3      = pars(40);
        h0      = pars(41);
        p0      = pars(42); %LdeP Total guess, can't find value in Gavina.

    % Drug parameters
    % Gavina
    E2exo           = pars(43);
    P4exo           = pars(44);


    % Calculate current aux equations
	E2=e0+e1*GrF+e2*DomF+e3*Lut4+E2exo;
	P4=p0+p1*Lut3+p2*Lut4+P4exo;
	Inh=h0+h1*DomF+h2*Lut2+h3*Lut3;

    % Calculate delayed aux equation
    ylag1 = lag_sol(:,1);
    IHdel=h0+h1*ylag1(7)+h2*ylag1(11)+h3*ylag1(12);


end
