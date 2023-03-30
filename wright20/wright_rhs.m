function dydt = wright_rhs(t, y, lag, pars)
 % Author: Andrew Wright, 2016
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
        % Pituitary parameters
      		v_0LH     = pars(1);
      		v_1LH     = pars(2);
      		Km_LH     = pars(3);
      		Ki_LH_P   = pars(4);
            k_LH      = pars(5);
            c_LH_P    = pars(6);
      		c_LH_E    = pars(7);
      		a_LH      = pars(8);
      		v_FSH     = pars(9);
            Ki_FSH_IH = pars(10);
            k_FSH     = pars(11);
            c_FSH_P   = pars(12);
            c_FSH_E   = pars(13);
      		a_FSH     = pars(14);    		
      		v         = pars(15);

        % Ovarian parameters
      		b         = pars(16);
      		Ki_RcF_P  = pars(17);
      		xi        = pars(18);
      		c_1       = pars(19);
      		c_2       = pars(20);
      		c_3       = pars(21);
      		c_4       = pars(22);
      		d_1       = pars(23);
      		d_2       = pars(24);
      		k_1       = pars(25);
      		k_2       = pars(26);
      		k_3       = pars(27);
      		k_4       = pars(28);
      		alpha     = pars(29);
      		gamma     = pars(30);

    % Get ovarian hormones
        [E2,P4,~,IHdel] = wright_OvarianHormones(y,lag,pars);

	% Synthesis and release of LH
        syn_LH=(v_0LH + v_1LH * Hill_plus(E2,Km_LH,8)) / (1+P4/Ki_LH_P);
		rel_LH=k_LH*RP_LH*((1+c_LH_P*P4)/(1+c_LH_E*E2));
		clear_LH=a_LH*LH;

	% Synthesis and release of FSH
		syn_FSH=v_FSH/(1+IHdel/Ki_FSH_IH);
		rel_FSH=k_FSH*(1+c_FSH_P*P4)*RP_FSH/(1+c_FSH_E*E2^2);
		clear_FSH=a_FSH*FSH;

	% Reserve pool and portal diff eq
		d_RP_FSH=syn_FSH-rel_FSH;
		d_FSH=1/v*rel_FSH-clear_FSH;

		d_RP_LH=syn_LH-rel_LH;
		d_LH=1/v*rel_LH-clear_LH;

	% Ovarian Stages
		d_RcF = (b*FSH+c_1*FSH*RcF)/(1+P4/Ki_RcF_P)^xi-c_2*(LH^alpha)*RcF;
		d_GrF = c_2*(LH^alpha)*RcF-c_3*LH*GrF;
		d_DomF = c_3*LH*GrF-c_4*LH^gamma*DomF;
		d_Sc1 = c_4*((LH)^gamma)*DomF-d_1*Sc1;
		d_Sc2 = d_1*Sc1-d_2*Sc2;
		d_Lut1 = d_2*Sc2-k_1*Lut1;
		d_Lut2 = k_1*Lut1-k_2*Lut2;
		d_Lut3 = k_2*Lut2-k_3*Lut3;
		d_Lut4 = k_3*Lut3-k_4*Lut4;

	% Assign values for return
		dydt = [d_RP_LH d_LH d_RP_FSH d_FSH d_RcF d_GrF d_DomF d_Sc1 d_Sc2 d_Lut1 d_Lut2 d_Lut3 d_Lut4]';

end
