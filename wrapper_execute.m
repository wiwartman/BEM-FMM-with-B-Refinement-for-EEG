clear all; %#ok<CLALL>
cd Model\
    model1_setup_base_model;
    model2_add_AMR;
cd ..\
bem1_setup_solution;
bem2_setup_integrals;
bem3_charge_engine;
bem3_surface_field_p;
bem3_surface_field_b;

