@everywhere function Control(t,x,data_dictionary)

  # Get the control parameters from the data dictionary -
  parameter_array = data_dictionary["control_parameter_array"]
  number_of_rates = data_dictionary["number_of_rates"]

  # Some of our species are parameters, they are encoded in experimental_input_array -
  experimental_input_array = data_dictionary["experimental_input_array"]
  TFPI              = experimental_input_array[1]
  FV                = experimental_input_array[2]
  FVIII             = experimental_input_array[3]
  FIX               = experimental_input_array[4]
  FX                = experimental_input_array[5]
  fXIII             = experimental_input_array[6]
  TAFI              = experimental_input_array[7]

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

  # Alias the control parameters -
  # Initiation
  alpha_trigger_activation        = parameter_array[1]
  order_trigger_activation        = parameter_array[2]
  alpha_trigger_inhibition_TFPI   = parameter_array[3]
  order_trigger_inhibition_TFPI   = parameter_array[4]

  # Amplification -
  alpha_amplification_APC         = parameter_array[5]
  order_amplification_APC         = parameter_array[6]
  alpha_amplification_TFPI        = parameter_array[7]
  order_amplification_TFPI        = parameter_array[8]

  # APC generation -
  alpha_shutdown_APC              = parameter_array[9]
  order_shutdown_APC              = parameter_array[10]

  alpha_fib_inh_fXIII             = parameter_array[11]
  order_fib_inh_fXIII             = parameter_array[12]
  alpha_fib_inh_TAFI              = parameter_array[13]
  order_fib_inh_TAFI              = parameter_array[14]

  # Control constants for plasmin activation term
  alpha_tPA_inh_PAI_1             = parameter_array[15]
  order_tPA_inh_PAI_1             = parameter_array[16]
  alpha_uPA_inh_PAI_1             = parameter_array[17]
  order_uPA_inh_PAI_1             = parameter_array[18]

  # Fomulate the control terms -
  control_array = ones(number_of_rates)

  # Initiation control -
  initiation_trigger_term   =   (alpha_trigger_activation*TRIGGER.^order_trigger_activation)/(1 + (alpha_trigger_activation*TRIGGER.^order_trigger_activation));
  initiation_TFPI_term      =   1 - (alpha_trigger_inhibition_TFPI*TFPI.^order_trigger_inhibition_TFPI)/(1 + (alpha_trigger_inhibition_TFPI*TFPI.^order_trigger_inhibition_TFPI));
  control_array[1]          =   minimum([initiation_trigger_term,initiation_TFPI_term])

  # Amplification FIIa feedback  -
  inhibition_term           =   1 - (alpha_amplification_APC*APC.^order_amplification_APC)/(1 + (alpha_amplification_APC*APC.^order_amplification_APC));
  inhibition_term_TFPI      =   1 - (alpha_amplification_TFPI*TFPI.^order_amplification_TFPI)/(1 +(alpha_amplification_TFPI*TFPI.^order_amplification_TFPI));
  factor_product            =   FV*FX*FVIII*FIX;
  factor_amplification_term =   (0.1*factor_product.^2)/(1+(0.1*factor_product.^2));
  control_array[2]          =   minimum([inhibition_term,inhibition_term_TFPI,factor_amplification_term])

  # Shutdown term -
  shutdown_term             =   (alpha_shutdown_APC*FIIa.^order_shutdown_APC)/(1 + (alpha_shutdown_APC*FIIa.^order_shutdown_APC));
  control_array[3]          =   shutdown_term

  # control_term_fibrinolysis_TAFI -
  control_term_fibrinolysis_TAFI  = 1 - ((alpha_fib_inh_TAFI.*TAFI).^order_fib_inh_TAFI)/(1+((alpha_fib_inh_TAFI.*TAFI).^order_fib_inh_TAFI))
  if (control_term_fibrinolysis_TAFI>0.001)
    control_array[10] = control_term_fibrinolysis_TAFI
  end

  control_term_fibrinolysis_TAFI  = 1 - ((alpha_fib_inh_TAFI.*TAFI).^order_fib_inh_TAFI)/(1+((alpha_fib_inh_TAFI.*TAFI).^order_fib_inh_TAFI))
  if (control_term_fibrinolysis_TAFI>0.001)
    control_array[11] = control_term_fibrinolysis_TAFI
  end

  # control_term_fibrinolysis_fXIII -
  control_term_fibrinolysis_fXIII  = 1 - ((alpha_fib_inh_fXIII.*fXIII).^order_fib_inh_fXIII)./(1+((alpha_fib_inh_fXIII.*fXIII).^order_fib_inh_fXIII))
  if (control_term_fibrinolysis_fXIII>0.001)
    control_array[12] = control_term_fibrinolysis_fXIII
  end

  # return -
  return control_array
end
