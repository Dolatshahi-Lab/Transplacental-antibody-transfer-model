function [diggdt] = diggdt(t,x,p)

% This model simulates the transport of IgG subclasses (1-4) through the
% placenta based on FcRn- or FcgRIIb-mediated transport as in STBs or ECs,
% respectively.
% Erdogan 2023; Dolatshahi Lab
% 02/23/2023

diggdt = zeros(35,1);

%assume maternal blood does not change (IgG regeneration)
diggdt(1) = 0; %free IgG in maternal blood
diggdt(2) = 0; %free IgG in maternal blood
diggdt(3) = 0; %free IgG in maternal blood
diggdt(4) = 0; %free IgG in maternal blood

diggdt(5) = (p.k_up*x(1) - p.k_deg*x(5) + p.koff1*x(9) - p.kon1*x(5)*x(13))/p.v_stb; %free IgG in STB endosome
diggdt(6) = (p.k_up*x(2) - p.k_deg*x(6) + p.koff2*x(10) - p.kon2*x(6)*x(13))/p.v_stb; %free IgG in STB endosomes
diggdt(7) = (p.k_up*x(3) - p.k_deg*x(7) + p.koff3*x(11) - p.kon3*x(7)*x(13))/p.v_stb; %free IgG in STB endosomes
diggdt(8) = (p.k_up*x(4) - p.k_deg*x(8) + p.koff4*x(12) - p.kon4*x(8)*x(13))/p.v_stb; %free IgG in STB endosomes

diggdt(9) = (-p.koff1*x(9) + p.kon1*x(5)*x(13) - p.k_t*x(9))/p.v_stb; %concentration of FcRn-IgG1 compley
diggdt(10) = (-p.koff2*x(10) + p.kon2*x(6)*x(13) - p.k_t*x(10))/p.v_stb; %concentration of FcRn-IgG2 compley
diggdt(11) = (-p.koff3*x(11) + p.kon3*x(7)*x(13) - p.k_t*x(11))/p.v_stb; %concentration of FcRn-IgG3 compley
diggdt(12) = (-p.koff4*x(12) + p.kon4*x(8)*x(13) - p.k_t*x(12))/p.v_stb; %concentration of FcRn-IgG4 compley

diggdt(13) =(-p.kon1*x(13)*x(5) + p.koff1*x(9) -p.kon2*x(13)*x(6) + p.koff2*x(10) ...
     - p.kon3*x(13)*x(7) + p.koff3*x(11) - p.kon4*x(13)*x(8) + p.koff4*x(12) ...
     + p.k_t*x(9) + p.k_t*x(10) + p.k_t*x(11) + p.k_t*x(12) + ...
     (2*p.fcrn_STB_curve.p1*t + p.fcrn_STB_curve.p2))/p.v_stb; %unbound FcRn in endosome

% %revised equations 8/21/2023 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diggdt(14) = (p.k_t*x(9) + p.koff1b*x(18) - p.kon1b*x(14)*x(22) - p.k_up*x(14))/p.v_str; %free IgG in stroma unbound IgG1
diggdt(15) = (p.k_t*x(10) + p.koff2b*x(19) - p.kon2b*x(15)*x(22) - p.k_up*x(15))/p.v_str; %free IgG in stroma
diggdt(16) = (p.k_t*x(11) + p.koff3b*x(20) - p.kon3b*x(16)*x(22) - p.k_up*x(16))/p.v_str; %free IgG in stroma
diggdt(17) = (p.k_t*x(12) + p.koff4b*x(21) - p.kon4b*x(17)*x(22) - p.k_up*x(17))/p.v_str; %free IgG in stroma

diggdt(18) = (-p.koff1b*x(18) + p.kon1b*x(14)*x(22) - p.k_t*x(18))/p.v_str; %concentration of FcgRIIb-IgG1 compley at the surface
diggdt(19) = (-p.koff2b*x(19) + p.kon2b*x(15)*x(22) - p.k_t*x(19))/p.v_str; %concentration of FcgRIIb-IgG2 compley at the surface
diggdt(20) = (-p.koff3b*x(20) + p.kon3b*x(16)*x(22) - p.k_t*x(20))/p.v_str; %concentration of FcgRIIb-IgG3 compley at the surface
diggdt(21) = (-p.koff4b*x(21) + p.kon4b*x(17)*x(22) - p.k_t*x(21))/p.v_str; %concentration of FcgRIIb-IgG4 compley at the surface

diggdt(31) = (-p.koff1*x(31) + p.kon1*x(23)*x(35) - p.k_t*x(31))/p.v_ec; %concentration of FcRn-IgG1 compley in endosomes
diggdt(32) = (-p.koff2*x(32) + p.kon2*x(24)*x(35) - p.k_t*x(32))/p.v_ec; %concentration of FcRn-IgG2 compley in endosomes
diggdt(33) = (-p.koff3*x(33) + p.kon3*x(25)*x(35) - p.k_t*x(33))/p.v_ec; %concentration of FcRn-IgG3 compley in endosomes
diggdt(34) = (-p.koff4*x(34) + p.kon4*x(26)*x(35) - p.k_t*x(34))/p.v_ec; %concentration of FcRn-IgG4 compley in endosomes

diggdt(22) = (-p.kon1b*x(22)*x(14) + p.koff1b*x(18) -p.kon2b*x(22)*x(15) + p.koff2b*x(19) ...
     - p.kon3b*x(22)*x(16) + p.koff3b*x(20) - p.kon4b*x(22)*x(17) + p.koff4b*x(21) ...
     + p.k_t*(x(18) + x(19) + x(20) + x(21)) + ...
     (2*p.fcgr2b_EC_curve.p1*t + p.fcgr2b_EC_curve.p2))/p.v_str; %unbound/free FcgRIIb on EC surface

diggdt(35) = (-p.kon1*x(35)*x(23) + p.koff1*x(31) -p.kon2*x(35)*x(24) + p.koff2*x(32) ...
     - p.kon3*x(35)*x(25) + p.koff3*x(33) - p.kon4*x(35)*x(26) + p.koff4*x(34) ...
     + p.k_t*(x(31) + x(32) + x(33) + x(34)) + ...
     (2*p.fcrn_EC_curve.p1*t + p.fcrn_EC_curve.p2))/p.v_ec; %unbound/free FcRn in EC endosomes

diggdt(23) = (p.k_up*x(14) - p.kon1*x(23)*x(35) + p.koff1*x(31) - p.k_deg*x(23))/p.v_ec; %free IgG in endothelial cells endosomes
diggdt(24) = (p.k_up*x(15) - p.kon2*x(24)*x(35) + p.koff2*x(32) - p.k_deg*x(24))/p.v_ec; %free IgG in endothelial cells endosomes
diggdt(25) = (p.k_up*x(16) - p.kon3*x(25)*x(35) + p.koff3*x(33) - p.k_deg*x(25))/p.v_ec; %free IgG in endothelial cells endosomes
diggdt(26) = (p.k_up*x(17) - p.kon4*x(26)*x(35) + p.koff4*x(34) - p.k_deg*x(26))/p.v_ec; %free IgG in endothelial cells endosomes

diggdt(27) = (p.k_t*(x(18)+x(31)) - p.d_ab*x(27))/p.v_f; %free IgG in fetal blood
diggdt(28) = (p.k_t*(x(19)+x(32)) - p.d_ab*x(28))/p.v_f; %free IgG in fetal blood
diggdt(29) = (p.k_t*(x(20)+x(33)) - p.d_ab*x(29))/p.v_f; %free IgG in fetal blood
diggdt(30) = (p.k_t*(x(21)+x(34)) - p.d_ab*x(30))/p.v_f; %free IgG in fetal blood

% %assume maternal blood does not change (IgG regeneration)
% diggdt(1) = 0; %maternal blood
% diggdt(2) = 0; %maternal blood
% diggdt(3) = 0; %maternal blood
% diggdt(4) = 0; %maternal blood
% 
% diggdt(5) = (p.k_up*x(1) - p.k_deg*x(5) + p.koff1*x(9) - p.kon1*x(5)*x(13))/p.v_stb; %STB endosomes unbound IgG1
% diggdt(6) = (p.k_up*x(2) - p.k_deg*x(6) + p.koff2*x(10) - p.kon2*x(6)*x(13))/p.v_stb; %STB endosomes
% diggdt(7) = (p.k_up*x(3) - p.k_deg*x(7) + p.koff3*x(11) - p.kon3*x(7)*x(13))/p.v_stb; %STB endosomes
% diggdt(8) = (p.k_up*x(4) - p.k_deg*x(8) + p.koff4*x(12) - p.kon4*x(8)*x(13))/p.v_stb; %STB endosomes
% 
% diggdt(9) = (-p.koff1*x(9) + p.kon1*x(5)*x(13) - p.k_t*x(9))/p.v_stb; %concentration of FcRn-IgG1 complex
% diggdt(10) = (-p.koff2*x(10) + p.kon2*x(6)*x(13)  - p.k_t*x(10))/p.v_stb; %concentration of FcRn-IgG2 complex
% diggdt(11) = (-p.koff3*x(11) + p.kon3*x(7)*x(13)  - p.k_t*x(11))/p.v_stb; %concentration of FcRn-IgG3 complex
% diggdt(12) = (-p.koff4*x(12) + p.kon4*x(8)*x(13)  - p.k_t*x(12))/p.v_stb; %concentration of FcRn-IgG4 complex
% 
% diggdt(13) =(-p.kon1*x(13)*x(5) + p.koff1*x(9) -p.kon2*x(13)*x(6) + p.koff2*x(10) ...
%      - p.kon3*x(13)*x(7) + p.koff3*x(11) - p.kon4*x(13)*x(8) + p.koff4*x(12) ...
%      + p.k_t*x(9) + p.k_t*x(10) + p.k_t*x(11) + p.k_t*x(12) + ...
%      (2*p.fcrn_EC_curve.p1*t + p.fcrn_EC_curve.p2))/p.v_stb; %unbound FcRn in endosome
% %revised equations 8/21/2023
% %checked
% diggdt(14) = (p.k_t*x(9) + p.koff1b*x(18) - p.kon1b*x(14)*x(22) - p.k_up*x(14))/p.v_str; %stroma unbound IgG1
% diggdt(15) = (p.k_t*x(10) + p.koff2b*x(19) - p.kon2b*x(15)*x(22) - p.k_up*x(15))/p.v_str; %stroma
% diggdt(16) = (p.k_t*x(11) + p.koff3b*x(20) - p.kon3b*x(16)*x(22) - p.k_up*x(16))/p.v_str; %stroma
% diggdt(17) = (p.k_t*x(12) + p.koff4b*x(21) - p.kon4b*x(17)*x(22) - p.k_up*x(17))/p.v_str; %stroma
% 
% %checked
% diggdt(18) = (-p.koff1b*x(18) + p.kon1b*x(14)*x(22) - p.k_t*x(18))/p.v_str; %concentration of FcgRIIb-IgG1 complex at the surface
% diggdt(19) = (-p.koff2b*x(19) + p.kon2b*x(15)*x(22) - p.k_t*x(19))/p.v_str; %concentration of FcgRIIb-IgG2 complex at the surface
% diggdt(20) = (-p.koff3b*x(20) + p.kon3b*x(16)*x(22) - p.k_t*x(20))/p.v_str; %concentration of FcgRIIb-IgG3 complex at the surface
% diggdt(21) = (-p.koff4b*x(21) + p.kon4b*x(17)*x(22) - p.k_t*x(21))/p.v_str; %concentration of FcgRIIb-IgG4 complex at the surface
% 
% %checked
% diggdt(31) = (-p.koff1*x(31) + p.kon1*x(23)*x(35) - p.k_t*x(31))/p.v_ec; %concentration of FcRn-IgG1 complex in endosomes
% diggdt(32) = (-p.koff2*x(32) + p.kon2*x(24)*x(35) - p.k_t*x(32))/p.v_ec; %concentration of FcRn-IgG2 complex in endosomes
% diggdt(33) = (-p.koff3*x(33) + p.kon3*x(25)*x(35) - p.k_t*x(33))/p.v_ec; %concentration of FcRn-IgG3 complex in endosomes
% diggdt(34) = (-p.koff4*x(34) + p.kon4*x(26)*x(35) - p.k_t*x(34))/p.v_ec; %concentration of FcRn-IgG4 complex in endosomes
% 
% diggdt(22) = (-p.kon1b*x(22)*x(14) + p.koff1b*x(18) -p.kon2b*x(22)*x(15) + p.koff2b*x(19) ...
%      - p.kon3b*x(22)*x(16) + p.koff3b*x(20) - p.kon4b*x(22)*x(17) + p.koff4b*x(21) ...
%      + p.k_t*(x(18) + x(19) + x(20) + x(21)) + ...
%      (2*p.fcgr2b_EC_curve.p1*t + p.fcgr2b_EC_curve.p2))/p.v_stb; %unbound/free FcgRIIb on EC surface
% 
% diggdt(35) = (-p.kon1*x(35)*x(23) + p.koff1*x(31) -p.kon2*x(35)*x(24) + p.koff2*x(32) ...
%      - p.kon3*x(35)*x(25) + p.koff3*x(33) - p.kon4*x(35)*x(26) + p.koff4*x(34) ...
%      + p.k_t*(x(31) + x(32) + x(33) + x(34)) + ...
%      (2*p.fcrn_EC_curve.p1*t + p.fcrn_EC_curve.p2))/p.v_ec; %unbound/free FcRn in EC endosomes
% 
% %checked
% diggdt(23) = (p.k_up*x(14) - p.kon1*x(23)*x(35) + p.koff1*x(31))/p.v_ec; %endothelial cells
% diggdt(24) = (p.k_up*x(15) - p.kon2*x(24)*x(35) + p.koff2*x(32))/p.v_ec; %endothelial cells
% diggdt(25) = (p.k_up*x(16) - p.kon3*x(25)*x(35) + p.koff3*x(33))/p.v_ec; %endothelial cells
% diggdt(26) = (p.k_up*x(17) - p.kon4*x(26)*x(35) + p.koff4*x(34))/p.v_ec; %endothelial cells
% 
% %checked
% diggdt(27) = (p.k_t*(x(18)+x(31)) - p.d_ab*x(27))/p.v_f; %fetal blood
% diggdt(28) = (p.k_t*(x(19)+x(32)) - p.d_ab*x(28))/p.v_f; %fetal blood
% diggdt(29) = (p.k_t*(x(20)+x(33)) - p.d_ab*x(29))/p.v_f; %fetal blood
% diggdt(30) = (p.k_t*(x(21)+x(34)) - p.d_ab*x(30))/p.v_f; %fetal blood

% original equations
% diggdt(14) = (p.k_t*x(9) + p.koff1b*x(18) - p.kon1b*x(14)*x(22) + p.koff1*x(31) - p.kon1*x(14)*x(35))/p.v_str; %stroma unbound IgG1
% diggdt(15) = (p.k_t*x(10) + p.koff2b*x(19) - p.kon2b*x(15)*x(22) + p.koff2*x(32) - p.kon2*x(15)*x(35))/p.v_str; %stroma
% diggdt(16) = (p.k_t*x(11) + p.koff3b*x(20) - p.kon3b*x(16)*x(22) + p.koff3*x(33) - p.kon3*x(16)*x(35))/p.v_str; %stroma
% diggdt(17) = (p.k_t*x(12) + p.koff4b*x(21) - p.kon4b*x(17)*x(22) + p.koff4*x(34) - p.kon4*x(17)*x(35))/p.v_str; %stroma
% 
% diggdt(18) = (-p.koff1b*x(18) + p.kon1b*x(14)*x(22) - p.k_up*x(18))/p.v_str; %concentration of FcgRIIb-IgG1 complex
% diggdt(19) = (-p.koff2b*x(19) + p.kon2b*x(15)*x(22)  - p.k_up*x(19))/p.v_str; %concentration of FcgRIIb-IgG2 complex
% diggdt(20) = (-p.koff3b*x(20) + p.kon3b*x(16)*x(22)  - p.k_up*x(20))/p.v_str; %concentration of FcgRIIb-IgG3 complex
% diggdt(21) = (-p.koff4b*x(21) + p.kon4b*x(17)*x(22)  - p.k_up*x(21))/p.v_str; %concentration of FcgRIIb-IgG4 complex
% 
% diggdt(31) = (-p.koff1*x(31) + p.kon1*x(14)*x(35) - p.k_up*x(31))/p.v_str; %concentration of FcRn-IgG1 complex
% diggdt(32) = (-p.koff2*x(32) + p.kon2*x(15)*x(35) - p.k_up*x(32))/p.v_str; %concentration of FcRn-IgG2 complex
% diggdt(33) = (-p.koff3*x(33) + p.kon3*x(16)*x(35) - p.k_up*x(33))/p.v_str; %concentration of FcRn-IgG3 complex
% diggdt(34) = (-p.koff4*x(34) + p.kon4*x(17)*x(35) - p.k_up*x(34))/p.v_str; %concentration of FcRn-IgG4 complex
% 
% diggdt(22) =(-p.kon1b*x(22)*x(14) + p.koff1b*x(18) -p.kon2b*x(22)*x(15) + p.koff2b*x(19) ...
%      - p.kon3b*x(22)*x(16) + p.koff3b*x(20) - p.kon4b*x(22)*x(17) + p.koff4b*x(21) ...
%      + p.k_up*x(18) + p.k_up*x(19) + p.k_up*x(20) + p.k_up*x(21) + ...
%      (2*p.fcgr2b_curve.p1*t + p.fcgr2b_curve.p2))/p.v_str; %unbound FcgRIIb on EC surface
% 
% diggdt(35) =(-p.kon1*x(35)*x(14) + p.koff1*x(31) -p.kon2*x(35)*x(15) + p.koff2*x(32) ...
%      - p.kon3*x(35)*x(16) + p.koff3*x(33) - p.kon4*x(35)*x(17) + p.koff4*x(34) ...
%      + p.k_up*x(31) + p.k_up*x(32) + p.k_up*x(33) + p.k_up*x(34) + ...
%      (2*p.fcrn_ec_curve.p1*t + p.fcrn_ec_curve.p2))/p.v_str; %unbound FcRn on EC surface
% 
% diggdt(23) = (p.k_up*(x(18)+x(31)) - p.k_t*x(23))/p.v_ec; %endothelial cells
% diggdt(24) = (p.k_up*(x(19)+x(32)) - p.k_t*x(24))/p.v_ec; %endothelial cells
% diggdt(25) = (p.k_up*(x(20)+x(33)) - p.k_t*x(25))/p.v_ec; %endothelial cells
% diggdt(26) = (p.k_up*(x(21)+x(34)) - p.k_t*x(26))/p.v_ec; %endothelial cells
% 
% diggdt(27) = (p.k_t*x(23) - p.d_ab*x(27))/p.v_f; %fetal blood
% diggdt(28) = (p.k_t*x(24) - p.d_ab*x(28))/p.v_f; %fetal blood
% diggdt(29) = (p.k_t*x(25) - p.d_ab*x(29))/p.v_f; %fetal blood
% diggdt(30) = (p.k_t*x(26) - p.d_ab*x(30))/p.v_f; %fetal blood

end