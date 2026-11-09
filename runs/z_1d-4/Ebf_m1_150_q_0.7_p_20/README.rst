Ebf_m1_150_q_0.7_p_20
=====================

150 + 105 solar masses, 20-day initial period, Z = 0.0001.
Pextra_factor = 1 (MESA default) with fixed atmosphere opacity.
No cool-wind enhancement (``x_ctrl(3) = 0``).
Start from the normal initial ZAMS models with ``./rn``.

Run
---

Use MESA branch ``EbF/improve_split_merge_amr_mesh`` at commit
``a4cb2640dffc9b3652eb7f40155eefd09dc561d3``. Initialize the MESA SDK and
set ``MESA_DIR`` to that checkout. Install MESA from its root::

    git lfs pull
    NPROCS=1 OMP_NUM_THREADS=1 ./install

Unpack this case into a new directory, then run there::

    make -j1
    OMP_NUM_THREADS=1 ./rn

Output
------

The wind calibration is stored separately for each star in native extras
state. Photos produced by this calculation preserve it for later restarts.
Use ``./re <photo>`` when restarting one of those new photos.

Profiles retain the original output cadence and include ``conv_vel``,
``mlt_vc`` (local MLT/TDC velocity), and ``effective_mixing_speed``
(``3*D_mix/ell`` from all mixing contributions), in cm/s.

Fresh initialization and two accepted binary steps were checked on one thread.

Notes
-----

- 2026-11-09: Copy of ``../Ebf_m1_150_q_0.7_p_20_Pextra1p5`` with
  ``Pextra_factor`` set to 1 and the cool-wind boost off, run from scratch
  with ``./rn``.
- 2026-11-09: MESA has been updated with the kap preprocessor bugfix (Type 1
  OPAL/OP tables returned spurious ``log10(kappa) = 0`` at low T / low rho when
  the interpolation stencil hit missing ``X = 0`` data). The fix is in
  ``kap/preprocessor/src/opal_type1.f`` (``if ((iadvance == 0) .and. (kmin <= 0)) return``
  and ``alt(max(1,kmin))``). The opacity tables were regenerated
  (``kap/kap_data.tar.xz``) and installed, and the kap cache was rebuilt. These are
  uncommitted local changes on top of commit ``a4cb2640d`` in ``$MESA_DIR``
  (see ``docs/source/known_bugs.rst``, ``opal_op_type1_preprocessor_bug``).
  This run uses the fixed tables; ``../Ebf_m1_150_q_0.7_p_20_Pextra1p5`` does not.
