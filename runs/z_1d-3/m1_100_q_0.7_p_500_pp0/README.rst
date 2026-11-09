.. _binary_ppisn_progenitor:

************************
binary_ppisn_progenitor
************************

.. tags:: star, binary, massive-star, mass-transfer, pair-instability

This development case constructs the binary-evolution progenitor used by
``dev_TDC_ppisn_from_binary``. It evolves an initially
:math:`100\,M_\odot + 70\,M_\odot` binary with an orbital period of 10 days
and metallicity :math:`Z = 0.0001`.

Evolution and model handoff
===========================

``inlist_binary`` selects ``inlist1`` for the initial donor and ``inlist2``
for the initial accretor. Both stars share the physics in ``inlist_both`` and
``inlist_extra``. Select each star's stopping condition with ``x_integer_ctrl(1)``
in ``inlist1`` and ``inlist2``:

* ``1``: central He-4 mass fraction at or below ``1d-8``.
* ``2`` (default): central ``log10(T/K)`` at or above ``9d0``.

Either star may finish first. At its first accepted model satisfying the
selected condition, ``run_binary_extras`` saves that star's configured final
model and ``final_profile.data``. The finished star becomes a point mass at
its current mass, and the other star continues to its own stopping condition.
Following Neev's handoff, the orbit is reset to a circular separation of
``100000 Rsun``, with Eddington-limited retention and radiation-corrected
transfer enabled. This is a prescribed handoff, not a supernova calculation.

Both stellar models remain available in restart photos, but only the
unfinished star is evolved. The run ends when no evolved star remains
unfinished, including runs started with a point-mass companion. Other
termination conditions, including the existing L2 checks, remain active.

The configured final model names are ``MASS1donor_final.mod`` and
``MASS2accretor_final.mod``. The latter is the default input model for
``dev_TDC_ppisn_from_binary``.

Run sequence
============

Running ``./rn`` selects ``inlist_binary`` explicitly with ``MESA_INLIST``.
That file also loads ``inlist_pgbinary``. Running ``./re`` restarts the most
recent complete set of binary and stellar photos with the same inlist.

The model uses TDC with turbulent energy in the energy equation,
``approx21_cr60_plus_co56.net``, the deBoer
:math:`{}^{12}\mathrm{C}(\alpha,\gamma){}^{16}\mathrm{O}` rate, and
``mass_fraction_limit_for_Skye = 1d-5``.

Last-Updated: 2026-08 by EbF

Notes
=====

- 2026-11-09: Rerun of ``m1_100_q_0.7_p_500_pp0`` (old output deleted),
  recreated from the current ``binary_ppisn_progenitor`` template with m1 = 100,
  m2 = 70, P = 500 d, Z = 1d-3, ``Pextra_factor = 1``, no rotation. All
  pp-chain reactions in ``approx21_cr60_plus_co56.net`` now have
  ``special_rate_factor = 0`` (``inlist_both``): ``rpp_to_he3``,
  ``r_he3_he3_to_h1_h1_he4``, ``r_he3_ag_be7``, ``r_be7_wk_li7``,
  ``r_be7_pg_b8``, ``r_h2_pg_he3``. The previous pp0 run only switched off
  ``rpp_to_he3`` and ``r_he3_he3_to_h1_h1_he4``, leaving the PP-II/III
  branch (``r_he3_ag_be7``) active. Run from scratch with ``./rn``.
- 2026-11-09: MESA has been updated with the kap preprocessor bugfix (Type 1
  OPAL/OP tables returned spurious ``log10(kappa) = 0`` at low T / low rho when
  the interpolation stencil hit missing ``X = 0`` data). The fix is in
  ``kap/preprocessor/src/opal_type1.f``; the opacity tables were regenerated
  and installed and the kap cache rebuilt. These are uncommitted local changes
  on top of commit ``a4cb2640d`` in ``$MESA_DIR`` (see
  ``docs/source/known_bugs.rst``, ``opal_op_type1_preprocessor_bug``).
  This run uses the fixed tables; the previous pp0 run did not.
