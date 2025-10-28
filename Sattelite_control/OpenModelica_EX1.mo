model Full_Model
  // Limit parameters
  parameter Real T_limMax = 300 "Body 1 maximum temperature";
  parameter Real T_limMin = 290 "Body 1 minimum temperature";
  parameter Real T_10perMax = 294.44415 "+10% temperature";
  parameter Real T_10perMin = 293.85585 "-10% temperature";
  // Control Parameters
  parameter Modelica.Units.SI.Resistance R = 0.1 "Motor resistance";
  parameter Modelica.Units.SI.Inductance L = 0.001 "Motor inductance";
  parameter Real k_m = 0.3 "Motor constant (Nm/A)";
  // Gain
  parameter Real k_p = 4e-5 "Gain constant (V/K)";
  // Radiator Mechanics
  parameter Modelica.Units.SI.Length L_r = 0.5  "Radiator length";
  parameter Modelica.Units.SI.Mass m_r = 0.2    "Radiator mass";
  parameter Modelica.Units.SI.Inertia J_r = (1/3)*m_r*L_r^2 "Radiator moment of inertia";
  // Radiator Emissivity
  parameter Modelica.Units.SI.Emissivity epsilon_min = 0.01 "Minimum emissivity of radiator";
  parameter Modelica.Units.SI.Emissivity epsilon_max = 0.98 "Maximum emissivity of radiator";
  // Radiator Areas
  parameter Modelica.Units.SI.Area A_r1 = 2*(1.5*0.5) + 2*(1.49*0.5) + (0.5*0.5)  "Body 1 radiator area";
  parameter Modelica.Units.SI.Area A_r2 = (0.95*0.5) + 2*(0.01*0.95) + (0.01*0.5) "Body 2 radiator area";
  parameter Modelica.Units.SI.Area A_r3 = (0.95*0.5) + 2*(0.01*0.95) + (0.01*0.5) "Body 3 radiator area";
  parameter Modelica.Units.SI.Area A_r4 = (0.5*0.5) "Body 4 radiator area";
  parameter Modelica.Units.SI.Area A_r5 = (0.5*0.5) "Body 5 radiator area";
  // Solar Panel Areas
  parameter Modelica.Units.SI.Area A_s1 = (0.5*0.5)   "Body 1 exposed to sun";
  parameter Modelica.Units.SI.Area A_s2 = (0.95*0.5)  "Body 2 (Solar panel)";
  parameter Modelica.Units.SI.Area A_s3 = (0.95*0.5)  "Body 3 (Solar panel)";
  // Thermal network parameters
  parameter Modelica.Units.SI.HeatCapacity C_1 = 1.5e5 "Main body heat capacity (J/K)";
  parameter Modelica.Units.SI.HeatCapacity C_2 = 1187.5 "Solar panel heat capacity (J/K)";
  parameter Modelica.Units.SI.HeatCapacity C_3 = 1187.5 "Solar panel heat capacity (J/K)";
  parameter Modelica.Units.SI.HeatCapacity C_4 = 30 "Radiator heat capacity (J/K)";
  parameter Modelica.Units.SI.HeatCapacity C_5 = 30 "Radiator heat capacity (J/K)";
  parameter Modelica.Units.SI.ThermalConductance G_12 = 10 "Conductance between node 1 and 2 (W/K)";
  parameter Modelica.Units.SI.ThermalConductance G_13 = 10 "Conductance between node 1 and 3 (W/K)";
  parameter Modelica.Units.SI.ThermalConductance G_14 = 10 "Conductance between node 1 and 4 (W/K)";
  parameter Modelica.Units.SI.ThermalConductance G_15 = 10 "Conductance between node 1 and 5 (W/K)";
  // Solar Radiation
  parameter Modelica.Units.SI.HeatFlux P_sun = 1350 "Solar radiation (W/m^2)";
  parameter Real alpha_1 = 0.6  "Absorptivity of body 1";
  parameter Real alpha_2 = 0.78 "Absorptivity of body 2 (solar panels)";
  parameter Real alpha_3 = 0.78 "Absorptivity of body 3 (solar panels)";
  parameter Modelica.Units.SI.Emissivity epsilon_1 = 0.45 "Emissivity of body 1";
  parameter Modelica.Units.SI.Emissivity epsilon_2 = 0.75 "Emissivity of body 2 (solar panels)";
  parameter Modelica.Units.SI.Emissivity epsilon_3 = 0.75 "Emissivity of body 3 (solar panels)";
  // Temperature and Angle Constraints
  parameter Modelica.Units.SI.Temperature T_ds(displayUnit = "K") = 3 "Deep Space temperature (K)";
  parameter Modelica.Units.SI.Temperature T_ref(displayUnit = "K") = 294.15 "Reference temperature (K)";
  parameter Modelica.Units.SI.Temperature T_start(displayUnit = "K") = 298.15 "Initial temperature (K)";
  parameter Modelica.Units.SI.Angle theta_min = -0.4*Modelica.Constants.pi "Radiator fully closed angle (rad)";
  parameter Modelica.Units.SI.Angle theta_max(displayUnit = "rad") = 0 "Radiator fully open angle (rad)";

  model emissivityVariation
  
  Modelica.Blocks.Interfaces.RealInput theta_input annotation(
      Placement(transformation(origin = {-134, 0}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-112, 0}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Interfaces.RealOutput eps_output annotation(
      Placement(transformation(origin = {126, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {126, 0}, extent = {{-10, -10}, {10, 10}})));
    parameter Real max_epsilon;
    parameter Real min_epsilon;
    constant Real pi = Modelica.Constants.pi;
  equation
eps_output = min_epsilon + ((max_epsilon - min_epsilon)/(0.4*pi))*(theta_input + 0.4*pi);
  annotation(
      Placement(transformation(origin = {4, -92}, extent = {{-10, -10}, {10, 10}})));
end emissivityVariation;

  model variableRadiation
  parameter Modelica.Units.SI.Area radiator_area;
    Modelica.Units.SI.HeatFlowRate rad_heat;
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(
      Placement(transformation(origin = {-114, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-114, 0}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b annotation(
      Placement(transformation(origin = {118, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {118, 0}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealInput epsilon_theta annotation(
      Placement(transformation(origin = {0, -120}, extent = {{-20, -20}, {20, 20}}, rotation = 90), iconTransformation(origin = {0, -104}, extent = {{-20, -20}, {20, 20}}, rotation = 90)));
  
  equation
  port_a.Q_flow = rad_heat;
  port_b.Q_flow = -rad_heat;
  rad_heat = epsilon_theta*radiator_area*Modelica.Constants.sigma*(port_a.T^4 - port_b.T^4);
  end variableRadiation;

  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor1(C = C_1, T(start = T_start, fixed = true))  annotation(
    Placement(transformation(origin = {172, 58}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor2(C = C_2, T(start = T_start, fixed = true))  annotation(
    Placement(transformation(origin = {104, 124}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor3(C = C_3, T(start = T_start, fixed = true))  annotation(
    Placement(transformation(origin = {238, 126}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor4(C = C_4, T(start = T_start, fixed = true))  annotation(
    Placement(transformation(origin = {108, -18}, extent = {{-10, 10}, {10, -10}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor5(C = C_5, T(start = T_start, fixed = true))  annotation(
    Placement(transformation(origin = {238, -20}, extent = {{-10, 10}, {10, -10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor12(G = G_12)  annotation(
    Placement(transformation(origin = {136, 96}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor13(G = G_13)  annotation(
    Placement(transformation(origin = {212, 98}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor14(G = G_14)  annotation(
    Placement(transformation(origin = {132, 14}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor15(G = G_15)  annotation(
    Placement(transformation(origin = {210, 14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation1(Gr = epsilon_1*A_r1)  annotation(
    Placement(transformation(origin = {160, -18}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation2(Gr = epsilon_2*A_r2)  annotation(
    Placement(transformation(origin = {66, 86}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation3(Gr = epsilon_3*A_r3)  annotation(
    Placement(transformation(origin = {278, 86}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow q_sun1(Q_flow = P_sun*alpha_1*A_s1)  annotation(
    Placement(transformation(origin = {186, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow q_sun2(Q_flow = P_sun*alpha_2*A_s2)  annotation(
    Placement(transformation(origin = {28, 114}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow q_sun3(Q_flow = P_sun*alpha_3*A_s3)  annotation(
    Placement(transformation(origin = {300, 116}, extent = {{10, -10}, {-10, 10}})));
  variableRadiation variableRadiation4(radiator_area = A_r4)  annotation(
    Placement(transformation(origin = {75, -51}, extent = {{-15, -15}, {15, 15}}, rotation = -90)));
  variableRadiation variableRadiation5(radiator_area = A_r5)  annotation(
    Placement(transformation(origin = {282, -50}, extent = {{-14, -14}, {14, 14}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature dsTemperature(T = T_ds)  annotation(
    Placement(transformation(origin = {176, -108}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  emissivityVariation emissivityVar(max_epsilon = epsilon_max, min_epsilon = epsilon_min)  annotation(
    Placement(transformation(origin = {-25, -99}, extent = {{-17, -17}, {17, 17}})));
  Modelica.Electrical.Analog.Basic.Ground ground annotation(
    Placement(transformation(origin = {-280, -82}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.Resistor resistor(R = R)  annotation(
    Placement(transformation(origin = {-236, -58}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.Inductor inductor(i(start = 0, fixed = true), L = L)  annotation(
    Placement(transformation(origin = {-198, -58}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.RotationalEMF emf(k = k_m)  annotation(
    Placement(transformation(origin = {-180, -96}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation(
    Placement(transformation(origin = {-250, -96}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Mechanics.Rotational.Components.Inertia radiatorInertia(J = J_r, phi(start = theta_min, fixed = true), w(start = 0, fixed = true), a(start = 0, fixed = true))  annotation(
    Placement(transformation(origin = {-146, -96}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Mechanics.Rotational.Sensors.AngleSensor angleSensor annotation(
    Placement(transformation(origin = {-118, -60}, extent = {{-10, -10}, {10, 10}})));

  model angleAndVoltageCondition
  // start paramters
  parameter Real condTheta_max;
  parameter Real condTheta_min;
  // end parameters
  Modelica.Blocks.Interfaces.RealInput condTheta_in annotation(
      Placement(transformation(origin = {-68, 118}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {-60, 116}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  Modelica.Blocks.Interfaces.RealInput condVoltage_in annotation(
      Placement(transformation(origin = {-108, -22}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-118, -44}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Interfaces.RealOutput condTheta_out annotation(
      Placement(transformation(origin = {110, 46}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {108, 50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput condVoltage_out annotation(
      Placement(transformation(origin = {110, -50}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, -44}, extent = {{-10, -10}, {10, 10}})));
  equation
  if condTheta_in > condTheta_max and condVoltage_in < 0 then
    condTheta_out = condTheta_max;
    condVoltage_out = 0;
  elseif condTheta_in < condTheta_min and condVoltage_in > 0 then
    condTheta_out = condTheta_min;
    condVoltage_out = 0;
  else
      condTheta_out = condTheta_in;
      condVoltage_out = condVoltage_in;
  end if;
  end angleAndVoltageCondition;
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor annotation(
    Placement(transformation(origin = {-212, 48}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Math.Feedback voltageFeedback annotation(
    Placement(transformation(origin = {-212, -126}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant refTemp(k = T_ref)  annotation(
    Placement(transformation(origin = {-212, -162}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
  Modelica.Blocks.Math.Gain kp_gain(k = -k_p)  annotation(
    Placement(transformation(origin = {-142, -142}, extent = {{-10, -10}, {10, 10}})));
  angleAndVoltageCondition condBlock(condTheta_max = theta_max, condTheta_min = theta_min)  annotation(
    Placement(transformation(origin = {-74, -106}, extent = {{-16, -16}, {16, 16}})));
equation
  connect(heatCapacitor2.port, thermalConductor12.port_b) annotation(
    Line(points = {{104, 114}, {104, 96}, {126, 96}}, color = {191, 0, 0}));
  connect(thermalConductor12.port_a, heatCapacitor1.port) annotation(
    Line(points = {{146, 96}, {152, 96}, {152, 48}, {172, 48}}, color = {191, 0, 0}));
  connect(thermalConductor13.port_a, heatCapacitor1.port) annotation(
    Line(points = {{202, 98}, {192, 98}, {192, 48}, {172, 48}}, color = {191, 0, 0}));
  connect(thermalConductor13.port_b, heatCapacitor3.port) annotation(
    Line(points = {{222, 98}, {238, 98}, {238, 116}}, color = {191, 0, 0}));
  connect(thermalConductor14.port_a, heatCapacitor1.port) annotation(
    Line(points = {{142, 14}, {152, 14}, {152, 48}, {172, 48}}, color = {191, 0, 0}));
  connect(thermalConductor15.port_a, heatCapacitor1.port) annotation(
    Line(points = {{200, 14}, {192, 14}, {192, 48}, {172, 48}}, color = {191, 0, 0}));
  connect(bodyRadiation1.port_a, heatCapacitor1.port) annotation(
    Line(points = {{160, -8}, {160, 28}, {172, 28}, {172, 48}}, color = {191, 0, 0}));
  connect(q_sun1.port, heatCapacitor1.port) annotation(
    Line(points = {{186, -8}, {186, 28}, {172, 28}, {172, 48}}, color = {191, 0, 0}));
  connect(bodyRadiation2.port_a, heatCapacitor2.port) annotation(
    Line(points = {{76, 86}, {104, 86}, {104, 114}}, color = {191, 0, 0}));
  connect(bodyRadiation3.port_a, heatCapacitor3.port) annotation(
    Line(points = {{268, 86}, {238, 86}, {238, 116}}, color = {191, 0, 0}));
  connect(q_sun2.port, heatCapacitor2.port) annotation(
    Line(points = {{38, 114}, {104, 114}}, color = {191, 0, 0}));
  connect(q_sun3.port, heatCapacitor3.port) annotation(
    Line(points = {{290, 116}, {238, 116}}, color = {191, 0, 0}));
  connect(thermalConductor14.port_b, heatCapacitor4.port) annotation(
    Line(points = {{122, 14}, {108, 14}, {108, -8}}, color = {191, 0, 0}));
  connect(thermalConductor15.port_b, heatCapacitor5.port) annotation(
    Line(points = {{220, 14}, {238, 14}, {238, -10}}, color = {191, 0, 0}));
  connect(variableRadiation4.port_a, heatCapacitor4.port) annotation(
    Line(points = {{75, -34}, {75, -8}, {108, -8}}, color = {191, 0, 0}));
  connect(heatCapacitor5.port, variableRadiation5.port_b) annotation(
    Line(points = {{238, -10}, {238, -11}, {282, -11}, {282, -33}}, color = {191, 0, 0}));
  connect(bodyRadiation1.port_b, dsTemperature.port) annotation(
    Line(points = {{160, -28}, {160, -44}, {176, -44}, {176, -98}}, color = {191, 0, 0}));
  connect(bodyRadiation2.port_b, dsTemperature.port) annotation(
    Line(points = {{56, 86}, {32, 86}, {32, -88}, {176, -88}, {176, -98}}, color = {191, 0, 0}));
  connect(variableRadiation4.port_b, dsTemperature.port) annotation(
    Line(points = {{75, -69}, {75, -88}, {176, -88}, {176, -98}}, color = {191, 0, 0}));
  connect(bodyRadiation3.port_b, dsTemperature.port) annotation(
    Line(points = {{288, 86}, {308, 86}, {308, -88}, {176, -88}, {176, -98}}, color = {191, 0, 0}));
  connect(variableRadiation5.port_a, dsTemperature.port) annotation(
    Line(points = {{282, -66}, {282, -88}, {176, -88}, {176, -98}}, color = {191, 0, 0}));
  connect(emissivityVar.eps_output, variableRadiation4.epsilon_theta) annotation(
    Line(points = {{-4, -99}, {16, -99}, {16, -51}, {59, -51}}, color = {0, 0, 127}));
  connect(emissivityVar.eps_output, variableRadiation5.epsilon_theta) annotation(
    Line(points = {{-4, -99}, {16, -99}, {16, -134}, {312, -134}, {312, -50}, {297, -50}}, color = {0, 0, 127}));
  connect(heatCapacitor1.port, temperatureSensor.port) annotation(
    Line(points = {{172, 48}, {-202, 48}}, color = {191, 0, 0}));
  connect(signalVoltage.n, resistor.p) annotation(
    Line(points = {{-250, -86}, {-250, -58}, {-246, -58}}, color = {0, 0, 255}));
  connect(resistor.n, inductor.p) annotation(
    Line(points = {{-226, -58}, {-208, -58}}, color = {0, 0, 255}));
  connect(inductor.n, emf.p) annotation(
    Line(points = {{-188, -58}, {-180, -58}, {-180, -86}}, color = {0, 0, 255}));
  connect(signalVoltage.p, emf.n) annotation(
    Line(points = {{-250, -106}, {-180, -106}}, color = {0, 0, 255}));
  connect(emf.flange, radiatorInertia.flange_a) annotation(
    Line(points = {{-170, -96}, {-156, -96}}));
  connect(radiatorInertia.flange_b, angleSensor.flange) annotation(
    Line(points = {{-136, -96}, {-132, -96}, {-132, -60}, {-128, -60}}));
  connect(ground.p, signalVoltage.n) annotation(
    Line(points = {{-280, -72}, {-250, -72}, {-250, -86}}, color = {0, 0, 255}));
  connect(temperatureSensor.T, voltageFeedback.u1) annotation(
    Line(points = {{-223, 48}, {-310, 48}, {-310, -126}, {-220, -126}}, color = {0, 0, 127}));
  connect(voltageFeedback.u2, refTemp.y) annotation(
    Line(points = {{-212, -134}, {-212, -151}}, color = {0, 0, 127}));
  connect(voltageFeedback.y, kp_gain.u) annotation(
    Line(points = {{-202, -126}, {-162, -126}, {-162, -142}, {-154, -142}}, color = {0, 0, 127}));
  connect(angleSensor.phi, condBlock.condTheta_in) annotation(
    Line(points = {{-106, -60}, {-84, -60}, {-84, -87}}, color = {0, 0, 127}));
  connect(kp_gain.y, condBlock.condVoltage_in) annotation(
    Line(points = {{-130, -142}, {-110, -142}, {-110, -113}, {-93, -113}}, color = {0, 0, 127}));
  connect(condBlock.condVoltage_out, signalVoltage.v) annotation(
    Line(points = {{-56, -113}, {-50, -113}, {-50, -194}, {-276, -194}, {-276, -96}, {-262, -96}}, color = {0, 0, 127}));
  connect(condBlock.condTheta_out, emissivityVar.theta_input) annotation(
    Line(points = {{-57, -98}, {-46.5, -98}, {-46.5, -99}, {-44, -99}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "4.0.0")),
  Diagram(coordinateSystem(extent = {{-320, 140}, {320, -200}}), graphics = {Rectangle(origin = {-169, -108}, extent = {{-121, 66}, {121, -66}}), Text(textColor = {0, 0, 255}, extent = {{-214, -36}, {-214, -36}}, textString = "Control Model"), Text(origin = {-168, -33}, textColor = {0, 0, 255}, extent = {{-48, 7}, {48, -7}}, textString = "Control Model"), Rectangle(origin = {163, 7}, extent = {{-153, 131}, {153, -131}}), Text(origin = {169, 148}, textColor = {170, 0, 0}, extent = {{-59, 8}, {59, -8}}, textString = "Thermal Model")}),
  version = "");
end Full_Model;
