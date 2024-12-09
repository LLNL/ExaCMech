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
