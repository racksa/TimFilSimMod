// seeding.hpp

// =============================================================================
// Include guard
#ifndef MY_SEEDING_FUNCTIONS_HEADER_INCLUDED
#define MY_SEEDING_FUNCTIONS_HEADER_INCLUDED

#include "config.hpp"

#if SURFACE_OF_REVOLUTION_BODIES

    void seed_blobs(double *const blob_references, double *const polar_dir_refs, double *const azi_dir_refs, double *const normal_refs);
    void seed_filaments(double *const filament_references, double *const polar_dir_refs, double *const azi_dir_refs, double *const normal_refs);

#endif

#endif // MY_SEEDING_FUNCTIONS_HEADER_INCLUDED
