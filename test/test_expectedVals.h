/**
 * @file test_expectedVals.h
 *
 * @brief Shared regression-value fragment (see `setup_base.h` for the general
 * `setup_*.h`/`test_*.h` inclusion pattern) `#include`d by `test_evptn.cxx` and
 * `test_updst.cxx`. Both drive an FCC single-crystal update to convergence and check
 * the result against these previously-recorded reference numbers -- so a diff here
 * would mean the underlying math/solver changed, not that these values are somehow
 * independently "correct" from first principles.
 *
 * Selects one of four sets of `expectedNFEvals`/`expectedGdotVal`/`expectedE2`/
 * `expectedQ1` values based on the including file's `KIN_TYPE` (and, for the default
 * `KIN_TYPE`, `XM_MUSHY`) build-time macro -- i.e. one reference solution per
 * slip-geometry/kinetics-model combination under test (`KIN_TYPE == 3`: BCC +
 * `Kin_KMBalD_TFF`; `== 2`: HCP + `Kin_HCP_A`; `== 1`: FCC + `Kin_KMBalD_FFF`;
 * default: FCC + `Kin_Voce`, further split by `XM_MUSHY`).
 *
 * @note `expectedNFEvals` is the raw solver function-evaluation count; `test_evptn.cxx`
 * checks it directly, while `test_updst.cxx` checks `expectedNFEvals + 1` since it also
 * requests the tangent stiffness matrix, which costs one extra evaluation.
 */
#if KIN_TYPE == 3

static const int   expectedNFEvals = 9;
static const double expectedGdotVal = -0.24401180817524;
static const double expectedE2 = 0.0023952767304669;
static const double expectedQ1 = 0.9996875162757;

#elif KIN_TYPE == 2

static const int   expectedNFEvals = 16;
static const double expectedGdotVal = 0.15704792600045;
static const double expectedE2 = 0.0067661223196391;
static const double expectedQ1 = 0.9996875162757;

#elif KIN_TYPE == 1

static const int   expectedNFEvals = 23;
static const double expectedGdotVal = -0.2398180885495;
static const double expectedE2 = 0.00407276458021;
static const double expectedQ1 = 0.999687516276;

#else

#if XM_MUSHY
static const int   expectedNFEvals = 13;
static const double expectedGdotVal = -0.4607285834886;
static const double expectedE2 = 0.00082808734102797;
static const double expectedQ1 = 0.9956165282270;
#else
static const int   expectedNFEvals = 18;
static const double expectedGdotVal = -0.2475346625929;
static const double expectedE2 = 0.0009861349707681;
static const double expectedQ1 = 0.999687516276;
#endif

#endif
