void filament::find_fitted_shape_s(){

      s_to_use = std::vector<Real>(NSEG);
      s_to_use[0] = 0.0;
      s_to_use[NSEG-1] = 1.0;


    // Bilinear interpolation
    Real phi_index = 0.5*phase/PI; // = phase/(2*pi)
    phi_index -= myfil_floor(phi_index); // Map this ratio into [0,1)
    phi_index *= (*s_to_use_ref_ptr).num_rows;

    int phi_index_int_lower_bound = int(phi_index); // Rounds towards 0.

    matrix svec;

    if (phi_index_int_lower_bound == (*s_to_use_ref_ptr).num_rows){

        // This can only be true if phi_index == s_to_use_ref.num_rows
        // I don't think we can ever actually enter here unless rounding error messes things up,
        // but better safe than sorry.
        svec = (*s_to_use_ref_ptr).get_row(0);

    } else {

        // Otherwise, we're safely in the interior
        const matrix lower = (*s_to_use_ref_ptr).get_row(phi_index_int_lower_bound);
        const matrix upper = (phi_index_int_lower_bound == (*s_to_use_ref_ptr).num_rows - 1) ? (*s_to_use_ref_ptr).get_row(0) : (*s_to_use_ref_ptr).get_row(phi_index_int_lower_bound+1);
        svec = lower + (upper - lower)*(phi_index - phi_index_int_lower_bound);

    }

    for (int n = 1; n < NSEG-1; n++){

        Real s_index = ((*s_to_use_ref_ptr).num_cols - 1)*n/Real(NSEG-1);
        int s_index_lower_bound = int(s_index);

        s_to_use[n] = svec(s_index_lower_bound) + (s_index - s_index_lower_bound)*(svec(s_index_lower_bound+1) - svec(s_index_lower_bound));

    }


    }