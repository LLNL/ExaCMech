/**
 * @file ECMech_slipgeom.h
 * @brief Umbrella header for the slip-geometry class family.
 *
 * Includes the common `SlipGeom` base template (ECMech_slipgeom_base.h) along with the
 * concrete crystal-structure implementations: FCC (ECMech_slipgeom_fcc.h), BCC
 * (ECMech_slipgeom_bcc.h), and HCP (ECMech_slipgeom_hcp.h). Include this header rather
 * than the individual `slipgeom/` headers when more than one crystal structure is
 * needed.
 */

#pragma once

#include "slipgeom/ECMech_slipgeom_base.h"
#include "slipgeom/ECMech_slipgeom_fcc.h"
#include "slipgeom/ECMech_slipgeom_bcc.h"
#include "slipgeom/ECMech_slipgeom_hcp.h"
