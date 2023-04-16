function [diggdt_vax] = diggdt_vax(t,x,p,vax_subclass_input)
%Author: Remziye Erdogan
%Date: 02/23/2023
%Here, there
%are two categories of IgG:  Total IgG subclasses, and vaccine-induced
%antigen-specific IgG subclasses. Here, we can ask questions about how the
%polarization of vaccine-induced IgG subclasses affects transfer
%efficiency.

diggdt_vax = zeros(63,1);

%assume maternal blood does not change (IgG regeneration)
diggdt_vax(1) = 0; %maternal blood
diggdt_vax(2) = 0; %maternal blood
diggdt_vax(3) = 0; %maternal blood
diggdt_vax(4) = 0; %maternal blood

diggdt_vax(5) = (p.k_up*x(1) - p.k_deg*x(5) + p.koff1*x(9) - p.kon1*x(5)*x(13))/p.v_se; %STB endosomes unbound IgG1
diggdt_vax(6) = (p.k_up*x(2) - p.k_deg*x(6) + p.koff2*x(10) - p.kon2*x(6)*x(13))/p.v_se; %STB endosomes
diggdt_vax(7) = (p.k_up*x(3) - p.k_deg*x(7) + p.koff3*x(11) - p.kon3*x(7)*x(13))/p.v_se; %STB endosomes
diggdt_vax(8) = (p.k_up*x(4) - p.k_deg*x(8) + p.koff4*x(12) - p.kon4*x(8)*x(13))/p.v_se; %STB endosomes

diggdt_vax(9) = (-p.koff1*x(9) + p.kon1*x(5)*x(13) - p.k_t*x(9))/p.v_se; %concentration of FcRn-IgG1 complex
diggdt_vax(10) = (-p.koff2*x(10) + p.kon2*x(6)*x(13)  - p.k_t*x(10))/p.v_se; %concentration of FcRn-IgG2 complex
diggdt_vax(11) = (-p.koff3*x(11) + p.kon3*x(7)*x(13)  - p.k_t*x(11))/p.v_se; %concentration of FcRn-IgG3 complex
diggdt_vax(12) = (-p.koff4*x(12) + p.kon4*x(8)*x(13)  - p.k_t*x(12))/p.v_se; %concentration of FcRn-IgG4 complex

diggdt_vax(13) =(-p.kon1*x(13)*x(5) + p.koff1*x(9) -p.kon2*x(13)*x(6) + p.koff2*x(10) ...
     - p.kon3*x(13)*x(7) + p.koff3*x(11) - p.kon4*x(13)*x(8) + p.koff4*x(12) ...
     + p.k_t*x(9) + p.k_t*x(10) + p.k_t*x(11) + p.k_t*x(12) + ...
     -p.kon1*x(13)*x(35) + p.koff1*x(39) -p.kon2*x(13)*x(36) + p.koff2*x(40) ...
     - p.kon3*x(13)*x(37) + p.koff3*x(41) - p.kon4*x(13)*x(38) + p.koff4*x(42) ...
     + p.k_t*x(39) + p.k_t*x(40) + p.k_t*x(41) + p.k_t*x(42) + ...
     (2*p.fcrn_curve.p1*t + p.fcrn_curve.p2))/p.v_se; %unbound FcRn in endosome

diggdt_vax(14) = (p.k_t*x(9) + p.koff1b*x(18) - p.kon1b*x(14)*x(22))/p.v_s; %stroma unbound IgG1
diggdt_vax(15) = (p.k_t*x(10) + p.koff2b*x(19) - p.kon2b*x(15)*x(22))/p.v_s; %stroma
diggdt_vax(16) = (p.k_t*x(11) + p.koff3b*x(20) - p.kon3b*x(16)*x(22))/p.v_s; %stroma
diggdt_vax(17) = (p.k_t*x(12) + p.koff4b*x(21) - p.kon4b*x(17)*x(22))/p.v_s; %stroma

diggdt_vax(18) = (-p.koff1b*x(18) + p.kon1b*x(14)*x(22) - p.k_up*x(18))/p.v_s; %concentration of FcgRIIb-IgG1 complex
diggdt_vax(19) = (-p.koff2b*x(19) + p.kon2b*x(15)*x(22)  - p.k_up*x(19))/p.v_s; %concentration of FcgRIIb-IgG2 complex
diggdt_vax(20) = (-p.koff3b*x(20) + p.kon3b*x(16)*x(22)  - p.k_up*x(20))/p.v_s; %concentration of FcgRIIb-IgG3 complex
diggdt_vax(21) = (-p.koff4b*x(21) + p.kon4b*x(17)*x(22)  - p.k_up*x(21))/p.v_s; %concentration of FcgRIIb-IgG4 complex

diggdt_vax(22) =(-p.kon1b*x(22)*x(14) + p.koff1b*x(18) -p.kon2b*x(22)*x(15) + p.koff2b*x(19) ...
     - p.kon3b*x(22)*x(16) + p.koff3b*x(20) - p.kon4b*x(22)*x(17) + p.koff4b*x(21) ...
     + p.k_up*x(18) + p.k_up*x(19) + p.k_up*x(20) + p.k_up*x(21) + ...
     (2*p.fcgr2b_curve.p1*t + p.fcgr2b_curve.p2) ...
     - p.kon1b*x(22)*x(44) + p.koff1b*x(48) -p.kon2b*x(22)*x(45) + p.koff2b*x(49) ...
     - p.kon3b*x(22)*x(46) + p.koff3b*x(50) - p.kon4b*x(22)*x(47) + p.koff4b*x(51) ...
     + p.k_up*x(48) + p.k_up*x(49) + p.k_up*x(50) + p.k_up*x(51))/p.v_s; %unbound FcgRIIb on EC surface

diggdt_vax(23) = (p.k_up*x(18) - p.k_t*x(23))/p.v_ee; %endothelial cells
diggdt_vax(24) = (p.k_up*x(19) - p.k_t*x(24))/p.v_ee; %endothelial cells
diggdt_vax(25) = (p.k_up*x(20) - p.k_t*x(25))/p.v_ee; %endothelial cells
diggdt_vax(26) = (p.k_up*x(21) - p.k_t*x(26))/p.v_ee; %endothelial cells

diggdt_vax(27) = (p.k_t*x(23) - p.d_ab*x(27))/p.v_f; %fetal blood
diggdt_vax(28) = (p.k_t*x(24) - p.d_ab*x(28))/p.v_f; %fetal blood
diggdt_vax(29) = (p.k_t*x(25) - p.d_ab*x(29))/p.v_f; %fetal blood
diggdt_vax(30) = (p.k_t*x(26) - p.d_ab*x(30))/p.v_f; %fetal blood

%%%%%%%%%%%%%%%%%Antigen-specific IgG transport%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if t <= p.tvax %before vaccine administration
        diggdt_vax(61) = 0; %antigen concentration
        diggdt_vax(62) = 0; %short-lived plasma cells
        diggdt_vax(63) = 0; %long-lived plasma cells
    else %after vaccine administration, introduce antigen
        diggdt_vax(61) = (-p.delta_ag*x(61))/p.v_m; %antigen concentration in maternal plasma
        diggdt_vax(62) = (p.rho*p.k_asc*x(61) - p.delta_s*x(62))/p.v_m; %short-lived plasma cells
        diggdt_vax(63) = ((1-p.rho)*p.k_asc*x(61) - p.delta_l*x(63))/p.v_m; %long-lived plasma cells
    end

    % This if statemnt determines which IgG subclass that vaccine induced.
    if double(vax_subclass_input) == 1
        diggdt_vax(31) = (p.k_ab*(x(62) + x(63)) - p.delta_ab*x(31))/p.v_m; %IgG1 concentration
        diggdt_vax(32) = 0; diggdt_vax(33) = 0; diggdt_vax(34) = 0; 
    elseif double(vax_subclass_input) == 2
        diggdt_vax(32) = (p.k_ab*(x(62) + x(63)) - p.delta_ab*x(32))/p.v_m; %IgG1 concentration
        diggdt_vax(31) = 0; diggdt_vax(33) = 0; diggdt_vax(34) = 0; 
    elseif double(vax_subclass_input) == 3
        diggdt_vax(33) = (p.k_ab*(x(62) + x(63)) - p.delta_ab*x(33))/p.v_m; %IgG1 concentration
        diggdt_vax(32) = 0; diggdt_vax(31) = 0; diggdt_vax(34) = 0; 
    elseif double(vax_subclass_input) == 4
        diggdt_vax(34) = (p.k_ab*(x(62) + x(63)) - p.delta_ab*x(34))/p.v_m; %IgG1 concentration
        diggdt_vax(32) = 0; diggdt_vax(33) = 0; diggdt_vax(31) = 0; 
    elseif strcmp(vax_subclass_input,'all proportional to bulk')
        diggdt_vax(31) = (0.6*p.k_ab*(x(62) + x(63)) - p.delta_ab*x(31))/p.v_m; %IgG1 concentration
        diggdt_vax(32) = (0.3*p.k_ab*(x(62) + x(63)) - p.delta_ab*x(32))/p.v_m; %IgG1 concentration
        diggdt_vax(33) = (0.05*p.k_ab*(x(62) + x(63)) - p.delta_ab*x(33))/p.v_m; %IgG1 concentration
        diggdt_vax(34) = (0.05*p.k_ab*(x(62) + x(63)) - p.delta_ab*x(34))/p.v_m; %IgG1 concentration
    elseif strcmp(vax_subclass_input,'all equal')
        diggdt_vax(31) = (0.25*p.k_ab*(x(62) + x(63)) - p.delta_ab*x(31))/p.v_m; %IgG1 concentration
        diggdt_vax(32) = (0.25*p.k_ab*(x(62) + x(63)) - p.delta_ab*x(32))/p.v_m; %IgG1 concentration
        diggdt_vax(33) = (0.25*p.k_ab*(x(62) + x(63)) - p.delta_ab*x(33))/p.v_m; %IgG1 concentration
        diggdt_vax(34) = (0.25*p.k_ab*(x(62) + x(63)) - p.delta_ab*x(34))/p.v_m; %IgG1 concentration
    end

diggdt_vax(35) = (p.k_up*x(31) - p.k_deg*x(35) + p.koff1*x(39) - p.kon1*x(35)*x(13))/p.v_se; %STB endosomes unbound IgG1
diggdt_vax(36) = (p.k_up*x(32) - p.k_deg*x(36) + p.koff2*x(40) - p.kon2*x(36)*x(13))/p.v_se; %STB endosomes
diggdt_vax(37) = (p.k_up*x(33) - p.k_deg*x(37) + p.koff3*x(41) - p.kon3*x(37)*x(13))/p.v_se; %STB endosomes
diggdt_vax(38) = (p.k_up*x(34) - p.k_deg*x(38) + p.koff4*x(42) - p.kon4*x(38)*x(13))/p.v_se; %STB endosomes

diggdt_vax(39) = (-p.koff1*x(39) + p.kon1*x(35)*x(13) - p.k_t*x(39))/p.v_se; %concentration of FcRn-IgG1 complex
diggdt_vax(40) = (-p.koff2*x(40) + p.kon2*x(36)*x(13)  - p.k_t*x(40))/p.v_se; %concentration of FcRn-IgG2 complex
diggdt_vax(41) = (-p.koff3*x(41) + p.kon3*x(37)*x(13)  - p.k_t*x(41))/p.v_se; %concentration of FcRn-IgG3 complex
diggdt_vax(42) = (-p.koff4*x(42) + p.kon4*x(38)*x(13)  - p.k_t*x(42))/p.v_se; %concentration of FcRn-IgG4 complex

diggdt_vax(44) = (p.k_t*x(39) + p.koff1b*x(48) - p.kon1b*x(44)*x(22))/p.v_s; %stroma unbound IgG1
diggdt_vax(45) = (p.k_t*x(40) + p.koff2b*x(49) - p.kon2b*x(45)*x(22))/p.v_s; %stroma
diggdt_vax(46) = (p.k_t*x(41) + p.koff3b*x(50) - p.kon3b*x(46)*x(22))/p.v_s; %stroma
diggdt_vax(47) = (p.k_t*x(42) + p.koff4b*x(51) - p.kon4b*x(47)*x(22))/p.v_s; %stroma

diggdt_vax(48) = (-p.koff1b*x(48) + p.kon1b*x(44)*x(22) - p.k_up*x(48))/p.v_s; %concentration of FcgRIIb-IgG1 complex
diggdt_vax(49) = (-p.koff2b*x(49) + p.kon2b*x(45)*x(22) - p.k_up*x(49))/p.v_s; %concentration of FcgRIIb-IgG2 complex
diggdt_vax(50) = (-p.koff3b*x(50) + p.kon3b*x(46)*x(22) - p.k_up*x(50))/p.v_s; %concentration of FcgRIIb-IgG3 complex
diggdt_vax(51) = (-p.koff4b*x(51) + p.kon4b*x(47)*x(22) - p.k_up*x(51))/p.v_s; %concentration of FcgRIIb-IgG4 complex

diggdt_vax(53) = (p.k_up*x(48) - p.k_t*x(53))/p.v_ee; %endothelial cells
diggdt_vax(54) = (p.k_up*x(49) - p.k_t*x(54))/p.v_ee; %endothelial cells
diggdt_vax(55) = (p.k_up*x(50) - p.k_t*x(55))/p.v_ee; %endothelial cells
diggdt_vax(56) = (p.k_up*x(51) - p.k_t*x(56))/p.v_ee; %endothelial cells

diggdt_vax(57) = (p.k_t*x(53) - p.d_ab*x(57))/p.v_f; %fetal blood
diggdt_vax(58) = (p.k_t*x(54) - p.d_ab*x(58))/p.v_f; %fetal blood
diggdt_vax(59) = (p.k_t*x(55) - p.d_ab*x(59))/p.v_f; %fetal blood
diggdt_vax(60) = (p.k_t*x(56) - p.d_ab*x(60))/p.v_f; %fetal blood

%Intracellular FcgRIIb binding%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diggdt_vax(44) = (p.k_t*x(39) - p.k_up*x(44))/p.v_s; %stroma unbound IgG1
% diggdt_vax(45) = (p.k_t*x(40) - p.k_up*x(45))/p.v_s; %stroma
% diggdt_vax(46) = (p.k_t*x(41) - p.k_up*x(46))/p.v_s; %stroma
% diggdt_vax(47) = (p.k_t*x(42) - p.k_up*x(47))/p.v_s; %stroma
% 
% diggdt_vax(48) = (-p.koff1b*x(48) + p.kon1b*x(53)*x(22) - p.k_t*x(48))/p.v_s; %concentration of FcgRIIb-IgG1 complex
% diggdt_vax(49) = (-p.koff2b*x(49) + p.kon2b*x(54)*x(22) - p.k_t*x(49))/p.v_s; %concentration of FcgRIIb-IgG2 complex
% diggdt_vax(50) = (-p.koff3b*x(50) + p.kon3b*x(55)*x(22) - p.k_t*x(50))/p.v_s; %concentration of FcgRIIb-IgG3 complex
% diggdt_vax(51) = (-p.koff4b*x(51) + p.kon4b*x(56)*x(22) - p.k_t*x(51))/p.v_s; %concentration of FcgRIIb-IgG4 complex
% 
% diggdt_vax(53) = (p.k_up*x(44) - p.kon1b*x(53)*x(22) + p.koff1b*x(48))/p.v_ee; %endothelial cells
% diggdt_vax(54) = (p.k_up*x(45) - p.kon2b*x(54)*x(22) + p.koff2b*x(49))/p.v_ee; %endothelial cells
% diggdt_vax(55) = (p.k_up*x(46) - p.kon3b*x(55)*x(22) + p.koff3b*x(50))/p.v_ee; %endothelial cells
% diggdt_vax(56) = (p.k_up*x(47) - p.kon4b*x(56)*x(22) + p.koff4b*x(51))/p.v_ee; %endothelial cells
% 
% diggdt_vax(57) = (p.k_t*x(48) - p.d_ab*x(57))/p.v_f; %fetal blood
% diggdt_vax(58) = (p.k_t*x(49) - p.d_ab*x(58))/p.v_f; %fetal blood
% diggdt_vax(59) = (p.k_t*x(50) - p.d_ab*x(59))/p.v_f; %fetal blood
% diggdt_vax(60) = (p.k_t*x(51) - p.d_ab*x(60))/p.v_f; %fetal blood

% diggdt_vax(22) =(-p.kon1b*x(22)*x(14) + p.koff1b*x(18) -p.kon2b*x(22)*x(15) + p.koff2b*x(19) ...
%      - p.kon3b*x(22)*x(16) + p.koff3b*x(20) - p.kon4b*x(22)*x(17) + p.koff4b*x(21) ...
%      + p.k_up*x(18) + p.k_up*x(19) + p.k_up*x(20) + p.k_up*x(21) + ...
%      (2*p.fcgr2b_curve.p1*t + p.fcgr2b_curve.p2) ...
%      - p.kon1b*x(22)*x(53) + p.koff1b*x(48) -p.kon2b*x(22)*x(54) + p.koff2b*x(49) ...
%      - p.kon3b*x(22)*x(55) + p.koff3b*x(50) - p.kon4b*x(22)*x(56) + p.koff4b*x(51) ...
%      + p.k_t*x(48) + p.k_t*x(49) + p.k_t*x(50) + p.k_t*x(51))/p.v_s; %unbound FcgRIIb on EC surface


end

