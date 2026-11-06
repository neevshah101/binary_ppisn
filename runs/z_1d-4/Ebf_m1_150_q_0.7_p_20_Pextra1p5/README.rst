Ebf_m1_150_q_0.7_p_20_Pextra1p5
===============================

150 + 105 solar masses, 20-day initial period, Z = 0.0001.
Pextra_factor = 1.5 with fixed atmosphere opacity.
Start from the normal initial ZAMS models with ``./rn``.

Run
---

Use MESA branch ``EbF/improve_split_merge_amr_mesh`` at commit
``5a177b45c259ab7f30ccde116e093184a697db5c``. Initialize the MESA SDK and
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
