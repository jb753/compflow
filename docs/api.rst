.. _api:

API reference
=============

This page defines the application--program interface for the compflow module.
The nomenclature section sets out the names and symbols for the quantites
calculated by the code, then the main functions are listed in three categories.
First, evaluation of a non-dimensional group given Mach number; second, the
reverse operation of finding the Mach number that corresponds to the value of
another quantity; and third, evaluation of derivatives with respect to Mach
number. This page concludes with a description of the lookup table function.

.. note::

    The library makes minimal attempts to validate input data, in the interests
    of speed. So, for example, feeding most functions a negative Mach number
    will silently produce bogus results. If this is not what you want, validate
    your own data before calling the library.

.. _nomen:

Nomenclature
************


The code deals in the following non-dimensional groups:

+---------------------------------+--------------------------------------------+---------------+
| Non-dimensional group           | Mathematical symbol                        | Variable name |
+=================================+============================================+===============+
| Mach number                     | :math:`\Ma`                                | `Ma`          |
+---------------------------------+--------------------------------------------+---------------+
| Ratio of specific heats         | :math:`\gamma = \dfrac{c_p}{c_v}`          | `ga`          |
+---------------------------------+--------------------------------------------+---------------+
| Stagnation temperature ratio    | :math:`\dfrac{T_0}{T}`                     | `To_T`        |
+---------------------------------+--------------------------------------------+---------------+
| Stagnation pressure ratio       | :math:`\dfrac{p_0}{p}`                     | `Po_P`        |
+---------------------------------+--------------------------------------------+---------------+
| Stagnation density ratio        | :math:`\dfrac{\rho_0}{\rho}`               | `rhoo_rho`    |
+---------------------------------+--------------------------------------------+---------------+
| Normalised velocity             | :math:`\dfrac{V}{\sqrt{c_pT_0}}`           | `V_cpTo`      |
+---------------------------------+--------------------------------------------+---------------+
| Impulse function                | :math:`\dfrac{F}{\dot{m}\sqrt{c_pT_0}}`    | `F_mcpTo`     |
+---------------------------------+--------------------------------------------+---------------+
| Normalised mass flow            | :math:`\dfrac{\dot{m}\sqrt{c_pT_0}}{Ap_0}` | `mcpTo_APo`   |
| (or flow function)              |                                            |               |
+---------------------------------+--------------------------------------------+---------------+
| Static normalised mass flow     | :math:`\dfrac{\dot{m}\sqrt{c_pT_0}}{Ap}`   | `mcpTo_AP`    |
+---------------------------------+--------------------------------------------+---------------+
| Area to choking area ratio      | :math:`\dfrac{A}{A_*}`                     | `A_Acrit`     |
+---------------------------------+--------------------------------------------+---------------+
| Post-shock Mach number          | :math:`\Ma_\mathrm{sh}`                    | `Mash`        |
+---------------------------------+--------------------------------------------+---------------+
| Shock stagnation pressure ratio | :math:`\dfrac{p_{0\mathrm{sh}}}{p_0}`      | `Posh_Po`     |
+---------------------------------+--------------------------------------------+---------------+

   
The following dimensional quantities are used in the definitions:

+-----------------+---------------------------------------------+
| Symbol          | Dimensional quantity                        |
+-----------------+---------------------------------------------+
| :math:`A`       | Flow area                                   |
+-----------------+---------------------------------------------+
| :math:`A_*`     | Choking flow area                           |
+-----------------+---------------------------------------------+
| :math:`c_p`     | Specific heat capacity at constant pressure |
+-----------------+---------------------------------------------+
| :math:`c_v`     | Specific heat capacity at constant volume   |
+-----------------+---------------------------------------------+
| :math:`F`       | Impulse                                     |
+-----------------+---------------------------------------------+
| :math:`\dot{m}` | Mass flow rate                              |
+-----------------+---------------------------------------------+
| :math:`p`       | Static pressure                             |
+-----------------+---------------------------------------------+
| :math:`p_0`     | Stagnation pressure                         |
+-----------------+---------------------------------------------+
| :math:`\rho`    | Static density                              |
+-----------------+---------------------------------------------+
| :math:`\rho_0`  | Stagnation density                          |
+-----------------+---------------------------------------------+
| :math:`T`       | Static temperature                          |
+-----------------+---------------------------------------------+
| :math:`T_0`     | Stagnation temperature                      |
+-----------------+---------------------------------------------+
| :math:`V`       | Velocity                                    |
+-----------------+---------------------------------------------+

Evaluation from Mach number
***************************

.. autosummary::

   compflow.To_T_from_Ma
   compflow.Po_P_from_Ma
   compflow.rhoo_rho_from_Ma
   compflow.V_cpTo_from_Ma
   compflow.mcpTo_APo_from_Ma
   compflow.mcpTo_AP_from_Ma
   compflow.F_mcpTo_from_Ma
   compflow.A_Acrit_from_Ma
   compflow.Mash_from_Ma
   compflow.Posh_Po_from_Ma

.. autofunction:: compflow.To_T_from_Ma
.. autofunction:: compflow.Po_P_from_Ma
.. autofunction:: compflow.rhoo_rho_from_Ma
.. autofunction:: compflow.V_cpTo_from_Ma
.. autofunction:: compflow.F_mcpTo_from_Ma
.. autofunction:: compflow.mcpTo_APo_from_Ma
.. autofunction:: compflow.mcpTo_AP_from_Ma
.. autofunction:: compflow.A_Acrit_from_Ma
.. autofunction:: compflow.Mash_from_Ma
.. autofunction:: compflow.Posh_Po_from_Ma


Inversion to Mach number
************************

.. autosummary::

   compflow.Ma_from_To_T
   compflow.Ma_from_Po_P
   compflow.Ma_from_rhoo_rho
   compflow.Ma_from_V_cpTo
   compflow.Ma_from_F_mcpTo
   compflow.Ma_from_mcpTo_APo
   compflow.Ma_from_mcpTo_AP
   compflow.Ma_from_A_Acrit
   compflow.Ma_from_Mash
   compflow.Ma_from_Posh_Po


.. autofunction:: compflow.Ma_from_To_T
.. autofunction:: compflow.Ma_from_Po_P
.. autofunction:: compflow.Ma_from_rhoo_rho
.. autofunction:: compflow.Ma_from_V_cpTo
.. autofunction:: compflow.Ma_from_F_mcpTo
.. autofunction:: compflow.Ma_from_mcpTo_APo
.. autofunction:: compflow.Ma_from_mcpTo_AP
.. autofunction:: compflow.Ma_from_A_Acrit
.. autofunction:: compflow.Ma_from_Mash
.. autofunction:: compflow.Ma_from_Posh_Po


Derivatives with respect to Mach number
***************************************

.. autosummary::

   compflow.der_To_T_from_Ma
   compflow.der_Po_P_from_Ma
   compflow.der_rhoo_rho_from_Ma
   compflow.der_V_cpTo_from_Ma
   compflow.der_F_mcpTo_from_Ma
   compflow.der_mcpTo_APo_from_Ma
   compflow.der_mcpTo_AP_from_Ma
   compflow.der_A_Acrit_from_Ma
   compflow.der_Mash_from_Ma
   compflow.der_Posh_Po_from_Ma

.. autofunction:: compflow.der_To_T_from_Ma
.. autofunction:: compflow.der_Po_P_from_Ma
.. autofunction:: compflow.der_rhoo_rho_from_Ma
.. autofunction:: compflow.der_V_cpTo_from_Ma
.. autofunction:: compflow.der_F_mcpTo_from_Ma
.. autofunction:: compflow.der_mcpTo_APo_from_Ma
.. autofunction:: compflow.der_mcpTo_AP_from_Ma
.. autofunction:: compflow.der_A_Acrit_from_Ma
.. autofunction:: compflow.der_Mash_from_Ma
.. autofunction:: compflow.der_Posh_Po_from_Ma

Lookup normalised mass flow
***************************

.. autofunction:: compflow.lookup_mcpTo_APo


Dimensional functions
*********************

.. autofunction:: compflow.static_from_stagnation
.. autofunction:: compflow.stagnation_from_static
.. autofunction:: compflow.change_frame
