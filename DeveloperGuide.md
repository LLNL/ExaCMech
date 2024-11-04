# Developer Guide
This is a working document and will likely need to be updated as further refactors occur in the future. Alternatively, if we realize things need better documentation then this will be updated.

# Name convention:
We currently don't have a consistent naming convention defined as we didn't create one beforehand. However, here is one we'll try and strive towards in future refactors.

Class and Struct types generally should follow a PascalCase convention aka `SuperCoolClassName`.

Variable names should generally be lower case with `_` between words. Although, we're still working on improving consistency in-regards to that convention across the library.

Class member variables, we generally want those to start `m_` and then `var_name` so they would be `m_var_name`. The `m_` here denotes a member variable. However, you will notice that for the `static constexpr int/size_t` variables we drop the `m_` and just have the name. The decision behind this was that these are static and constant properties of the class, and don't represent parameters that the user is able to define at runtime or that might be derived from runtime defined parameters.

Struct member variables differ from class member variables as they are generally just data containers without any functions. So, we don't require the names to look like `m_var_name` but we instead expect them to be `var_name`.

Function names should also be lowercase with `_` between names. However as earlier noted, we still have not exactly enforced this at this point.

# Variable Names:

We have a wide range of different variables used through out the code. This is a short guide on at least some of them and is by no means a complete guide to everything in the library.

Next, we try to add the following substrings to names if the variables refer to quantities that either are some form of a deviatoric, deviatoric + pressure, or some variation of those vector quantities in Voigt notation. `d5` refers to a deviatoric vector that makes use of the actual 5D space representation. `d6` refers to a deviatoric vector that has values in the expanded 6D space. If the `d5` or `d6` has a character following that it refers to either the pressure (`p`) or volume (`v`), and these vectors are concatenation of the deviatoric and volumetric/pressure terms. So, a variable with a `d6p` substring would be of length 7 where the first 6 terms are the deviatoric terms and the final is the pressure term.

First, `elast` suffix names refer to the elastic strain and generally this involves quantities in the crystal/lattice frame.

`def_rate` is the deforamtion rate which is the symmetric portion of the velocity gradient.

`spin_vec` is `vect(W)` where W is the skew portion of the velocity gradient.

`rel_vol_ratios` is maybe not the best named array but it's made up of what can be thought of terms related to the relative volume of the point your dealing with where the 2 first indices (0 & 1) correspond to the old and new relative volumes respectively. Here, a relative volume > 1 signals material under volumetric expansion and a relative volume < 1 signals material under volumetric compression. Traditionally, one might set this to begin as a value of 1 but that is not a hard requirement.

Terms with `quat` in the name generally refer to the unit quaternion that maps the rotation of the starting crystal lattice reference frame to the current sample frame.

Terms that look like `h_state` refer to the hardening state.

Terms that have the name `gdot` refer to the macroscopic shearing rate per slip system of the material commonly seen written as $\dot{\gamma}$ in the literature.

`tkelv` is the temperature which is most commonly provided in Kelvin.

Terms with the suffix `_n` commonly refer to the begining of time step information.

Terms with the suffix `_u` or `_f` commonly refer to the end of time step information.

# Further Notes:
If you find something unclear in either the documentation or code please open an issue asking for clarification. Alternatively, if you find issues with those same topics feel free to open a PR with a fix as generally our developer team's bandwidth is limited so some of these items will be lower on the priority list for us.

