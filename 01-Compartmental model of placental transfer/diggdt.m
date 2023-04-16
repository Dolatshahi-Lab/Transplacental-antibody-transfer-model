function [diggdt] = diggdt(t,x,p)

% This model simulates the transport of IgG subclasses (1-4) through the
% placenta based on FcRn- or FcgRIIb-mediated transport as in STBs or ECs,
% respectively.
% Erdogan 2023; Dolatshahi Lab
% 02/23/2023

diggdt = zeros(30,1);

% diggdt(1) = -p.k_up*x(1)/p.v_m; %maternal blood
% diggdt(2) = -p.k_up*x(2)/p.v_m; %maternal blood
% diggdt(3) = -p.k_up*x(3)/p.v_m; %maternal blood
% diggdt(4) = -p.k_up*x(4)/p.v_m; %maternal blood
%assume maternal blood does not change (IgG regeneration)
diggdt(1) = 0; %maternal blood
diggdt(2) = 0; %maternal blood
diggdt(3) = 0; %maternal blood
diggdt(4) = 0; %maternal blood

diggdt(5) = (p.k_up*x(1) - p.k_deg*x(5) + p.koff1*x(9) - p.kon1*x(5)*x(13))/p.v_se; %STB endosomes unbound IgG1
diggdt(6) = (p.k_up*x(2) - p.k_deg*x(6) + p.koff2*x(10) - p.kon2*x(6)*x(13))/p.v_se; %STB endosomes
diggdt(7) = (p.k_up*x(3) - p.k_deg*x(7) + p.koff3*x(11) - p.kon3*x(7)*x(13))/p.v_se; %STB endosomes
diggdt(8) = (p.k_up*x(4) - p.k_deg*x(8) + p.koff4*x(12) - p.kon4*x(8)*x(13))/p.v_se; %STB endosomes

diggdt(9) = (-p.koff1*x(9) + p.kon1*x(5)*x(13) - p.k_t*x(9))/p.v_se; %concentration of FcRn-IgG1 complex
diggdt(10) = (-p.koff2*x(10) + p.kon2*x(6)*x(13)  - p.k_t*x(10))/p.v_se; %concentration of FcRn-IgG2 complex
diggdt(11) = (-p.koff3*x(11) + p.kon3*x(7)*x(13)  - p.k_t*x(11))/p.v_se; %concentration of FcRn-IgG3 complex
diggdt(12) = (-p.koff4*x(12) + p.kon4*x(8)*x(13)  - p.k_t*x(12))/p.v_se; %concentration of FcRn-IgG4 complex

diggdt(13) =(-p.kon1*x(13)*x(5) + p.koff1*x(9) -p.kon2*x(13)*x(6) + p.koff2*x(10) ...
     - p.kon3*x(13)*x(7) + p.koff3*x(11) - p.kon4*x(13)*x(8) + p.koff4*x(12) ...
     + p.k_t*x(9) + p.k_t*x(10) + p.k_t*x(11) + p.k_t*x(12) + ...
     (2*p.fcrn_curve.p1*t + p.fcrn_curve.p2))/p.v_se; %unbound FcRn in endosome

diggdt(14) = (p.k_t*x(9) + p.koff1b*x(18) - p.kon1b*x(14)*x(22))/p.v_s; %stroma unbound IgG1
diggdt(15) = (p.k_t*x(10) + p.koff2b*x(19) - p.kon2b*x(15)*x(22))/p.v_s; %stroma
diggdt(16) = (p.k_t*x(11) + p.koff3b*x(20) - p.kon3b*x(16)*x(22))/p.v_s; %stroma
diggdt(17) = (p.k_t*x(12) + p.koff4b*x(21) - p.kon4b*x(17)*x(22))/p.v_s; %stroma

diggdt(18) = (-p.koff1b*x(18) + p.kon1b*x(14)*x(22) - p.k_up*x(18))/p.v_s; %concentration of FcgRIIb-IgG1 complex
diggdt(19) = (-p.koff2b*x(19) + p.kon2b*x(15)*x(22)  - p.k_up*x(19))/p.v_s; %concentration of FcgRIIb-IgG2 complex
diggdt(20) = (-p.koff3b*x(20) + p.kon3b*x(16)*x(22)  - p.k_up*x(20))/p.v_s; %concentration of FcgRIIb-IgG3 complex
diggdt(21) = (-p.koff4b*x(21) + p.kon4b*x(17)*x(22)  - p.k_up*x(21))/p.v_s; %concentration of FcgRIIb-IgG4 complex

diggdt(22) =(-p.kon1b*x(22)*x(14) + p.koff1b*x(18) -p.kon2b*x(22)*x(15) + p.koff2b*x(19) ...
     - p.kon3b*x(22)*x(16) + p.koff3b*x(20) - p.kon4b*x(22)*x(17) + p.koff4b*x(21) ...
     + p.k_up*x(18) + p.k_up*x(19) + p.k_up*x(20) + p.k_up*x(21) + ...
     (2*p.fcgr2b_curve.p1*t + p.fcgr2b_curve.p2))/p.v_s; %unbound FcgRIIb on EC surface

diggdt(23) = (p.k_up*x(18) - p.k_t*x(23))/p.v_ee; %endothelial cells
diggdt(24) = (p.k_up*x(19) - p.k_t*x(24))/p.v_ee; %endothelial cells
diggdt(25) = (p.k_up*x(20) - p.k_t*x(25))/p.v_ee; %endothelial cells
diggdt(26) = (p.k_up*x(21) - p.k_t*x(26))/p.v_ee; %endothelial cells

diggdt(27) = (p.k_t*x(23) - p.d_ab*x(27))/p.v_f; %fetal blood
diggdt(28) = (p.k_t*x(24) - p.d_ab*x(28))/p.v_f; %fetal blood
diggdt(29) = (p.k_t*x(25) - p.d_ab*x(29))/p.v_f; %fetal blood
diggdt(30) = (p.k_t*x(26) - p.d_ab*x(30))/p.v_f; %fetal blood

end