"""This module wraps the fortran backend."""

import compflow_fort_from_Ma as fort_from_Ma
import compflow_fort_der_from_Ma as fort_der_from_Ma
import compflow_fort_to_Ma as fort_to_Ma

# Functions from Ma
To_T_from_Ma = fort_from_Ma.to_t
Po_P_from_Ma = fort_from_Ma.po_p
rhoo_rho_from_Ma = fort_from_Ma.rhoo_rho
V_cpTo_from_Ma = fort_from_Ma.v_cpto
mcpTo_APo_from_Ma = fort_from_Ma.mcpto_apo
mcpTo_AP_from_Ma = fort_from_Ma.mcpto_ap
A_Acrit_from_Ma = fort_from_Ma.a_acrit
Mash_from_Ma = fort_from_Ma.mash
Posh_Po_from_Ma = fort_from_Ma.posh_po

# Derivatives from Ma
der_To_T_from_Ma = fort_der_from_Ma.to_t
der_Po_P_from_Ma = fort_der_from_Ma.po_p
der_rhoo_rho_from_Ma = fort_der_from_Ma.rhoo_rho
der_V_cpTo_from_Ma = fort_der_from_Ma.v_cpto
der_mcpTo_APo_from_Ma = fort_der_from_Ma.mcpto_apo
der_mcpTo_AP_from_Ma = fort_der_from_Ma.mcpto_ap
der_A_Acrit_from_Ma = fort_der_from_Ma.a_acrit
der_Mash_from_Ma = fort_der_from_Ma.mash
der_Posh_Po_from_Ma = fort_der_from_Ma.posh_po

# Inversions to Ma
Ma_from_To_T = fort_to_Ma.to_t
Ma_from_Po_P = fort_to_Ma.po_p
Ma_from_rhoo_rho = fort_to_Ma.rhoo_rho
Ma_from_V_cpTo = fort_to_Ma.v_cpto
Ma_from_mcpTo_APo = fort_to_Ma.mcpto_apo
Ma_from_mcpTo_AP = fort_to_Ma.mcpto_ap
Ma_from_A_Acrit = fort_to_Ma.a_acrit
Ma_from_Mash = fort_to_Ma.mash
Ma_from_Posh_Po = fort_to_Ma.posh_po
