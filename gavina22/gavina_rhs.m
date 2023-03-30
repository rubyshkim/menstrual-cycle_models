function dydt = gavina_rhs(t, y, lag, pars)
 %LdeP Updated for Gavina 2022 December 29, 2022
 % Author: Gavina, de los Reyes, Olufsen, Lenhart, Ottesen, 2022
 % Purpose: model for use with dde23 solver, contains the differential equations        

 % % Input
 % 	- time
 % 	- current solutions
 % 	- lag solution vector
 % 	- parameter vector
 % % Returns
 % 	- column vector of equation values for use with solver

 	% get current state variable values
		RP_LH  = y(1); %current RP_LH level
		LH     = y(2); %current LH level
		RP_FSH = y(3);
		FSH    = y(4);
		RcF    = y(5); %current MsF level
		GrF    = y(6); %current GrF level
		DomF   = y(7); %current DomF level
		Sc1    = y(8); %current Sc1 level
		Sc2    = y(9); %current Sc2 level
		Lut1   = y(10); %current Lut1 level
		Lut2   = y(11); %current Lut2 level
		Lut3   = y(12); %current Lut3 level
		Lut4   = y(13); %current Lut4 level

	% Get parameters
        k_LH    = pars(1);
        alpha_LH= pars(2);
        V_0LH   = pars(3);
        V_1LH   = pars(4);
        km_LH   = pars(5);
        %----- Constant Inhibition Parameters
        ki_LHP  = pars(6);
        c_LHE   = pars(7);
        c_FSHE  = pars(8);
        K_iFSHInh=pars(9);
        wstar   = pars(10);
        qstar   = pars(11);
        %----- Constant for Stimulation Parameters
        c_LHP   = pars(12);
        c_FSHP  = pars(13);
        %----- 
        V_FSH   = pars(14);
        alpha_FSH=pars(15);
        k_FSH   = pars(16);
        tau     = pars(17);
        v       = pars(18);
        %----- Follicles transition factor from
        b       = pars(19);
        c2      = pars(20);
        c3      = pars(21);
        c4      = pars(22);
        d1      = pars(23);
        d2      = pars(24);
        k1      = pars(25);
        k2      = pars(26);
        k3      = pars(27);
        c1      = pars(28);
        k4      = pars(29);
        %-----
        alpha   = pars(30);
        gamma   = pars(31);
        e0      = pars(32);
        %----- Contribution factor parameters
        e1      = pars(33);
        e2      = pars(34);
        e3      = pars(35);
        p1      = pars(36);
        p2      = pars(37);
        h1      = pars(38);
        h2      = pars(39);
        h3      = pars(40);
        %-----
        h0      = pars(41);
        p0      = pars(42); %LdeP Total guess, can't find value in Gavina.
        %----- Treatments
        E2exo   = pars(43);
        P4exo   = pars(44);
    %%%%%%%%%%%%%%%%%%%%%%%

    % Get ovarian hormones
        [E2,P4,~,IHdel] = gavina_OvarianHormones(y,lag,pars);

        %%%%%%%%LdeP Stopped HERE%%%%%%%%%%

	% Synthesis and release of LH
        syn_LH=(V_0LH + V_1LH * Hill_plus(E2,km_LH,8)) / (1+P4/ki_LHP);
		rel_LH=k_LH*RP_LH*((1+c_LHP*P4)/(1+c_LHE*E2));
		clear_LH=alpha_LH*LH;

	% Synthesis and release of FSH
		syn_FSH=V_FSH/(1+IHdel/K_iFSHInh + P4/wstar); %LdeP P4/wstar differs from Wright
		rel_FSH=k_FSH*(1+c_FSHP*P4)*RP_FSH/(1+c_FSHE*E2^2);
		clear_FSH=alpha_FSH*FSH;

	% Reserve pool and portal diff eq
		d_RP_FSH=syn_FSH-rel_FSH;
		d_FSH=rel_FSH*(1/v)-clear_FSH;

		d_RP_LH=syn_LH-rel_LH;
		d_LH=rel_LH*(1/v)-clear_LH;

	% Ovarian Stages
		d_RcF = (b*FSH+c1*FSH*RcF)/(1+P4/qstar)-c2*(LH^alpha)*RcF; %qstar not in Wright
		d_GrF = c2*(LH^alpha)*RcF-c3*LH*GrF;
		d_DomF = c3*LH*GrF-c4*LH^gamma*DomF;
		d_Sc1 = c4*((LH)^gamma)*DomF-d1*Sc1;
		d_Sc2 = d1*Sc1-d2*Sc2;
		d_Lut1 = d2*Sc2-k1*Lut1;
		d_Lut2 = k1*Lut1-k2*Lut2;
		d_Lut3 = k2*Lut2-k3*Lut3;
		d_Lut4 = k3*Lut3-k4*Lut4;

	% Assign values for return
		dydt = [d_RP_LH d_LH d_RP_FSH d_FSH d_RcF d_GrF d_DomF d_Sc1 d_Sc2 d_Lut1 d_Lut2 d_Lut3 d_Lut4]';

end
