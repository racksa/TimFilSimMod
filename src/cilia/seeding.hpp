// seeding.hpp

// =============================================================================
// Include guard
#ifndef MY_SEEDING_FUNCTIONS_HEADER_INCLUDED
#define MY_SEEDING_FUNCTIONS_HEADER_INCLUDED

#include "../../config.hpp"

#if SURFACE_OF_REVOLUTION_BODIES or ROD or RIGIDWALL

    void seed_blobs(Real *const blob_references, Real *const polar_dir_refs, Real *const azi_dir_refs, Real *const normal_refs);
    void seed_filaments(Real *const filament_references, Real *const polar_dir_refs, Real *const azi_dir_refs, Real *const normal_refs);

    // void seed_rod_blobs(Real *const blob_references, Real *const polar_dir_refs, Real *const azi_dir_refs, Real *const normal_refs);

    void check_seeding(Real *const filament_references, Real *const polar_dir_refs, Real *const azi_dir_refs, Real *const normal_refs);


#endif

#endif // MY_SEEDING_FUNCTIONS_HEADER_INCLUDED
