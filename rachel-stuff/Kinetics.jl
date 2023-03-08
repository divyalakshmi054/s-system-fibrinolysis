@everywhere function Kinetics(t,x,data_dictionary)

  # Get the parameters from the data dictionary -
  parameter_array = data_dictionary["kinetic_parameter_array"]
  number_of_rates = data_dictionary["number_of_rates"]

  delay_term = 1 - 1/(1+exp(0.1*(t-60)))
  #delay_term = 1.0

  # Alias the species -
  FII                = x[1]
  FIIa               = x[2]
  PC                 = x[3]
  APC                = x[4]
  ATIII              = x[5]
  TM                 = x[6]
  TRIGGER            = x[7]

  Fibrin             = x[8]
  Plasmin            = x[9]
  Fibrinogen         = x[10]
  Plasminogen        = x[11]
  tPA                = x[12]
  uPA                = x[13]
  Fibrin_monomer     = x[14]
  Protofibril        = x[15]
  antiplasmin        = x[16]
  PAI_1              = x[17]
  Fiber              = x[18]

  # Some of our species are parameters, they are encoded in experimental_input_array -
  experimental_input_array = data_dictionary["experimental_input_array"]
  TFPI              = experimental_input_array[1]
  FV                = experimental_input_array[2]
  FVIII             = experimental_input_array[3]
  FIX               = experimental_input_array[4]
  FX                = experimental_input_array[5]
  fXIII             = experimental_input_array[6]
  TAFI              = experimental_input_array[7]

  # Alias the parameters -
  # Kinetic parameters for thrombin generation and regulation
  k_trigger                                 = parameter_array[1]
  K_FII_trigger                             = parameter_array[2]
  k_amplification                           = parameter_array[3]
  K_FII_amplification                       = parameter_array[4]
  k_APC_formation                           = parameter_array[5]
  K_PC_formation                            = parameter_array[6]
  k_inhibition_ATIII                        = parameter_array[7]

  # Kinetic parameters for fibrin generation
  k_cat_Fibrinogen                          = parameter_array[8]
  Km_Fibrinogen                             = parameter_array[9]
  k_fibrin_monomer_association              = parameter_array[10]
  k_protofibril_association                 = parameter_array[11]
  k_protofibril_monomer_association         = parameter_array[12]

  # Kinetic parameters for fibrinolysis
  k_cat_plasminogen_Fibrin_tPA              = parameter_array[13]
  Km_plasminogen_Fibrin_tPA                 = parameter_array[14]
  k_cat_plasminogen_Fibrin_uPA              = parameter_array[15]
  Km_plasminogen_Fibrin_uPA                 = parameter_array[16]
  k_cat_Fibrin                              = parameter_array[17]
  Km_Fibrin                                 = parameter_array[18]
  k_inhibit_plasmin                         = parameter_array[19]
  k_inhibit_PAI_1_APC                       = parameter_array[20]
  k_inhibit_tPA_PAI_1                       = parameter_array[21]
  k_inhibit_uPA_PAI_1                       = parameter_array[22]

  # Kinetic parameter for fibrin formation
  k_fibrin_formation                        = parameter_array[23]

  # Fiber degradation
  k_cat_fiber                               = parameter_array[24]
  Km_fiber                                  = parameter_array[25]

  # Plasmin activation
  Ka_plasminogen_Fibrin_tPA                 = parameter_array[26]
  K_plasminogen_Fibrin_tPA                  = parameter_array[27]
  k_cat_fibrinogen_deg                      = parameter_array[28]
  Km_fibrinogen_deg                         = parameter_array[29]

  F= Fibrin + Fiber;
  k_cat_plasminogen_tPA = k_cat_plasminogen_Fibrin_tPA*(F)/(K_plasminogen_Fibrin_tPA+F);
  Km_plasminogen_tPA = Km_plasminogen_Fibrin_tPA*(Ka_plasminogen_Fibrin_tPA+F)/(K_plasminogen_Fibrin_tPA+F);

  # k_cat_plasminogen_tPA=k_cat_plasminogen_Fibrin_tPA
  # Km_plasminogen_tPA=Km_plasminogen_Fibrin_tPA

  # Formulate the kinetic rates -
  rate_array      =   zeros(number_of_rates)
  rate_array[1]   =   (k_trigger*TRIGGER*(FII/(K_FII_trigger + FII)))*delay_term
  rate_array[2]   =   k_amplification*FIIa*(FII/(K_FII_amplification + FII));
  rate_array[3]   =   k_APC_formation*TM*(PC/(K_PC_formation + PC));
  rate_array[4]   =   k_inhibition_ATIII*(ATIII)*(FIIa^1.26);
  rate_array[5]   =   (FIIa*k_cat_Fibrinogen*Fibrinogen)/(Km_Fibrinogen+Fibrinogen)                             # Cleavage of fibrinopeptides A and/or B to form fibrin monomer
  rate_array[6]   =   k_fibrin_monomer_association*(Fibrin_monomer^2);                                          # Protofibril formation through association of fibrin monomers
  rate_array[7]   =   k_protofibril_association*(Protofibril^2);                                                # Association of protofibril-protofibril to form fibers
  rate_array[8]   =   k_protofibril_monomer_association*(Protofibril)*(Fibrin_monomer);                         # Fibril association with monomer forms protofibril again
  rate_array[9]   =   k_fibrin_formation*(Fiber);                                                               # Fibrin growth
  rate_array[10]  =   (tPA*k_cat_plasminogen_tPA*Plasminogen)/(Km_plasminogen_tPA+Plasminogen);
  rate_array[11]  =   (uPA*k_cat_plasminogen_Fibrin_uPA*Plasminogen)/(Km_plasminogen_Fibrin_uPA+Plasminogen);
  rate_array[12]  =   (Plasmin*k_cat_Fibrin*Fibrin)/(Km_Fibrin+Fibrin)
  rate_array[13]  =   k_inhibit_plasmin*Plasmin*(antiplasmin);
  rate_array[14]  =   k_inhibit_PAI_1_APC*(APC)*(PAI_1);
  rate_array[15]  =   k_inhibit_tPA_PAI_1*(tPA)*PAI_1;
  rate_array[16]  =   k_inhibit_uPA_PAI_1*(uPA)*PAI_1;
  rate_array[17]  =   (Plasmin*k_cat_fiber*Fiber)/(Km_fiber+Fiber)                                           # Plasmin degrading fiber
  rate_array[18]  =   (Plasmin*k_cat_fibrinogen_deg*Fibrinogen)/(Km_fibrinogen_deg+Fibrinogen)

  # return the rate array -
  return rate_array
end
