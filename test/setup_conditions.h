double def_rate_d6v_sample[ecmech::nsvp] = { -0.5, -0.5, 1.0,
                                      0.0, 0.0, 0.0,
                                      0.0 };
vecsVsa<ecmech::nsvp>(def_rate_d6v_sample, ecmech::sqr2b3); // scale so that def_rate_d5_sample comes out to unit magnitude
//
double def_rate_d5_sample[ecmech::ntvec];
svecToVecd(def_rate_d5_sample, def_rate_d6v_sample);

double dt = 1e-1;

double spin_vec_sample[ecmech::nwvec] = { 0.0, 0.0, 0.5 };

double rel_vol_ratios[ecmech::nvr] = { 1.0, 1.0, 0.0, 0.0 };
