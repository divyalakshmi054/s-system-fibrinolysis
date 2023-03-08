@everywhere function DataFile()

  # Initialize the data dictionary -
  data_dictionary = Dict()

  # setup system dimension -
  number_of_rates = 18
  number_of_species = 18

  # set experimental inputs ============================================================== #
  experimental_input_array = zeros(7)
  experimental_input_array[1] = 2.5           # 1 TFPI
  experimental_input_array[2] = 20            # 2 FV
  experimental_input_array[3] = 0.7           # 3 FVIII
  experimental_input_array[4] = 90.0          # 4 FIX
  experimental_input_array[5] = 170           # 5 FX
  experimental_input_array[6] = 93            # 6 FXIII
  experimental_input_array[7] = 70            # 7 TAFI
  # ====================================================================================== #

  # set the initial condition ============================================================ #
  initial_condition_array = zeros(number_of_species)
  initial_condition_array[1]  = 1400          # FII                = x[1]
  initial_condition_array[2]  = 0.0           # FIIa               = x[2]
  initial_condition_array[3]  = 60.0          # PC                 = x[3]
  initial_condition_array[4]  = 0.0           # APC                = x[4]
  initial_condition_array[5]  = 3400          # ATIII              = x[5]
  initial_condition_array[6]  = 1.0           # TM                 = x[6]
  initial_condition_array[7]  = 0.0           # TRIGGER            = x[7]
  initial_condition_array[8]  = 0.0           # Fibrin             = x[8]
  initial_condition_array[9]  = 0.0           # Plasmin            = x[9]
  initial_condition_array[10] = 10000         # Fibrinogen         = x[10]
  initial_condition_array[11] = 2000          # Plasminogen        = x[11]
  initial_condition_array[12] = 8.0           # tPA                = x[12]
  initial_condition_array[13] = 0.0           # uPA                = x[13]
  initial_condition_array[14] = 0.0           # Fibrin_monomer     = x[14]
  initial_condition_array[15] = 0.0           # Protofibril        = x[15]
  initial_condition_array[16] = 1180          # antiplasmin        = x[16]
  initial_condition_array[17] = 0.56          # PAI_1              = x[17]
  initial_condition_array[18] = 0.0           # Fiber              = x[18]
  # ===================================================================================== #

  # set the kinetic rate constants  ===================================================== #
  kinetic_parameter_array = Float64[]
  push!(kinetic_parameter_array,130.0)            # 1 k_trigger
  push!(kinetic_parameter_array,1400.0)           # 2 K_FII_trigger
  push!(kinetic_parameter_array,0.0001)           # 3 k_amplification
  push!(kinetic_parameter_array,100.0)            # 4 K_FII_amplification
  push!(kinetic_parameter_array,0.001)            # 5 k_APC_formation
  push!(kinetic_parameter_array,60.0)             # 6 K_PC_formation
  push!(kinetic_parameter_array,0.0000001)         # 7 k_inhibition_ATIII

  # Kinetic parameters for fibrin generation
  push!(kinetic_parameter_array,0.01)            # 8 k_cat_Fibrinogen
  push!(kinetic_parameter_array,1000.0)         # 9 Km_Fibrinogen
  push!(kinetic_parameter_array,0.01)         # 10 k_fibrin_monomer_association
  push!(kinetic_parameter_array,0.01)         # 11 k_protofibril_association
  push!(kinetic_parameter_array,0.01)         # 12 k_protofibril_monomer_association

  # Kinetic parameters for fibrinolysis
  push!(kinetic_parameter_array,2.0)           # 13 k_cat_plasminogen_Fibrin_tPA
  push!(kinetic_parameter_array,10.0)          # 14 Km_plasminogen_Fibrin_tPA
  push!(kinetic_parameter_array,0.1)           # 15 k_cat_plasminogen_Fibrin_uPA
  push!(kinetic_parameter_array,10.0)          # 16 Km_plasminogen_Fibrin_uPA
  push!(kinetic_parameter_array,1.0)           # 17 k_cat_Fibrin
  push!(kinetic_parameter_array,10.0)          # 18 Km_Fibrin
  push!(kinetic_parameter_array,0.01)          # 19 k_inhibit_plasmin
  push!(kinetic_parameter_array,0.01)          # 20 k_inhibit_PAI_1_APC
  push!(kinetic_parameter_array,0.01)          # 21 k_inhibit_tPA_PAI_1
  push!(kinetic_parameter_array,0.01)          # 22 k_inhibit_uPA_PAI_1

  # Kinetic parameter for fibrin formation
  push!(kinetic_parameter_array,0.01)          # 23 k_fibrin_formation

  # Fiber degradation
  push!(kinetic_parameter_array,100.0)          # 24 k_cat_fiber
  push!(kinetic_parameter_array,10.0)          # 25 Km_fiber

  # Plasmin activation
  push!(kinetic_parameter_array,0.01)          # 26 Ka_plasminogen_Fibrin_tPA
  push!(kinetic_parameter_array,0.01)          # 27 K_plasminogen_Fibrin_tPA
  push!(kinetic_parameter_array,1.00)          # 28 k_cat_fibrinogen_deg
  push!(kinetic_parameter_array,10.0)          # 29 Km_fibrinogen_deg
  # ==================================================================================== #

  # set the control rate constants ===================================================== #
  control_parameter_array = Float64[]
  push!(control_parameter_array,0.1)          # 1 alpha_trigger_activation
  push!(control_parameter_array,2.0)          # 2 order_trigger_activation
  push!(control_parameter_array,0.1)          # 3 alpha_trigger_inhibition_TFPI
  push!(control_parameter_array,2.0)          # 4 order_trigger_inhibition_TFPI

  # Amplification -
  push!(control_parameter_array,0.1)          # 5 alpha_amplification_APC
  push!(control_parameter_array,2.0)          # 6 order_amplification_APC
  push!(control_parameter_array,0.1)          # 7 alpha_amplification_TFPI
  push!(control_parameter_array,2.0)          # 8 order_amplification_TFPI

  # APC generation -
  push!(control_parameter_array,0.1)          # 9 alpha_shutdown_APC
  push!(control_parameter_array,2.0)          # 10 order_shutdown_APC

  # Control constants for TAFI inhibition -
  push!(control_parameter_array,0.1)          # 11 alpha_fib_inh_fXIII
  push!(control_parameter_array,2.0)          # 12 order_fib_inh_fXIII
  push!(control_parameter_array,0.1)          # 13 alpha_fib_inh_TAFI
  push!(control_parameter_array,2.0)          # 14 order_fib_inh_TAFI

  # Control constants for plasmin activation term
  push!(control_parameter_array,0.1)          # 15 alpha_tPA_inh_PAI_1
  push!(control_parameter_array,2.0)          # 16 order_tPA_inh_PAI_1
  push!(control_parameter_array,0.1)          # 17 alpha_uPA_inh_PAI_1
  push!(control_parameter_array,2.0)          # 18 order_uPA_inh_PAI_1
  # ==================================================================================== #

  # ===== DO NOT EDIT (BELOW) ============================================= #
  data_dictionary["initial_condition_array"] = initial_condition_array
  data_dictionary["number_of_rates"] = number_of_rates
  data_dictionary["kinetic_parameter_array"] = kinetic_parameter_array
  data_dictionary["control_parameter_array"] = control_parameter_array
  data_dictionary["experimental_input_array"] = experimental_input_array
  return data_dictionary
  # ===== DO NOT EDIT (ABOVE) ============================================= #

end
