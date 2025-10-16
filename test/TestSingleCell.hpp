#include <cxxtest/TestSuite.h>
#include "AbstractCvodeCell.hpp"
#include "CellProperties.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RegularStimulus.hpp"
#include "HeartConfig.hpp"
//#include "Shannon2004Cvode.hpp"
//#include "carro_2011_epiCvodeOpt.hpp"
/*
#include "carro_2011_epiCvodeOpt.hpp"
#include "carro_2011_epi_pig_fkatp_v3CvodeOpt.hpp"
#include "carro_2011_endo_pig_fkatp_v3CvodeOpt.hpp"
#include "carro_2011_m_pig_fkatp_v3CvodeOpt.hpp"

#include "carro_2011_epi_pig_fkatp_v4CvodeOpt.hpp"
#include "carro_2011_endo_pig_fkatp_v4CvodeOpt.hpp"
#include "carro_2011_m_pig_fkatp_v4CvodeOpt.hpp"

#include "carro_2011_epi_pig_fkatp_v4_nalCvodeOpt.hpp"

#include "carro_2011_epi_pig_fkatp_v4_nal_ito2gCvodeOpt.hpp"
#include "carro_2011_endo_pig_fkatp_v4_nal_ito2gCvodeOpt.hpp"
#include "carro_2011_m_pig_fkatp_v4_nal_ito2gCvodeOpt.hpp"
*/
//#include "carro_2011_epi_pig_fkatp_v4_nal_ito2g2caclmod_k1mod_ssCvodeOpt.hpp"
//#include "carro_2011_endo_pig_fkatp_v4_nal_ito2g2caclmod_k1mod_ssCvodeOpt.hpp"
//#include "carro_2011_m_pig_fkatp_v4_nal_ito2g2caclmod_k1mod_ssCvodeOpt.hpp"

//#include "carro_2011_epi_pig_fkatp_v4_nal_ito2g2caclmod_k1mod_ncCvodeOpt.hpp"
//#include "carro_2011_endo_pig_fkatp_v4_nal_ito2g2caclmod_k1mod_ncCvodeOpt.hpp"

//#include "carro_2011_epi_pig_fkatp_v4_nal_ito2g2caclmod_k1mod_nc_mod2CvodeOpt.hpp"
//#include "carro_2011_endo_pig_fkatp_v4_nal_ito2g2caclmod_k1mod_nc_mod2CvodeOpt.hpp"

//#include "carro_2011_epi_pig_scratch_3_go1_ssCvodeOpt.hpp"
//#include "carro_2011_endo_pig_scratch_3_go1_ssCvodeOpt.hpp"
//#include "carro_2011_m_pig_scratch_3_go1_ssCvodeOpt.hpp"

#include "carro_2011_epi_pig_scratch_4_go1CvodeOpt.hpp"
#include "carro_2011_endo_pig_scratch_4_go1CvodeOpt.hpp"
#include "carro_2011_m_pig_scratch_4_go1CvodeOpt.hpp"

//#include "pigexx_chasteCvodeOpt.hpp"

#include "SteadyStateRunner.hpp"
#include "FakePetscSetup.hpp"

class TestSingleCellSimulationTutorial : public CxxTest::TestSuite
{
public:
    void TestShannonSimulation() /*throw (Exception)*/
    {
#ifdef CHASTE_CVODE

        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(-50,1.0,1000.0,10.0)); //mag was 50.. threshhol -20
        //boost::shared_ptr<RegularStimulus> p_stimulus;
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        //boost::shared_ptr<AbstractCvodeCell> p_model(new Cellcarro_2011_epi_pig_fkatp_v3FromCellMLCvodeOpt(p_solver, p_stimulus));

        //boost::shared_ptr<AbstractCvodeCell> p_model(new Cellcarro_2011_epi_pig_fkatp_v4FromCellMLCvodeOpt(p_solver, p_stimulus));
        //boost::shared_ptr<AbstractCvodeCell> p_model(new Cellcarro_2011_endo_pig_fkatp_v4FromCellMLCvodeOpt(p_solver, p_stimulus));
        //boost::shared_ptr<AbstractCvodeCell> p_model(new Cellcarro_2011_m_pig_fkatp_v4FromCellMLCvodeOpt(p_solver, p_stimulus));

        //boost::shared_ptr<AbstractCvodeCell> p_model(new Cellcarro_2011_m_pig_fkatp_v4_nal_ito2gFromCellMLCvodeOpt(p_solver, p_stimulus));
        //boost::shared_ptr<AbstractCvodeCell> p_model(new Cellpigexx_chasteFromCellMLCvodeOpt(p_solver, p_stimulus));

        boost::shared_ptr<AbstractCvodeCell> p_model(new Cellcarro_2011_epi_pig_scratch_4_go1FromCellMLCvodeOpt(p_solver, p_stimulus));
        // adp90 225 -epi, 235 endo: gks 0.045 for endo, 255 m: + gkr 0.055
        
        // stim starts at 10 (was 0), dur = 1, ampl = 40, perdiod = 1000
        //boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();

        //p_regular_stim->SetPeriod(1000); //500 it effects the properties 

        std::string layer = "epi";

        p_model->SetMaxSteps(1e6);

        double abs_t = 1e-8;
        double rel_t = 1e-8;
        p_model->SetTolerances(abs_t, rel_t); // e-8 good

        // control
        
        //p_model->SetParameter("membrane_L_type_calcium_current_conductance", 1.0); //1.2
        //p_model->SetParameter("membrane_fast_sodium_current_conductance", 19.0); //20 // 36
        //p_model->SetParameter("extracellular_potassium_concentration", 5.0); // 5.9
        //p_model->SetParameter("extracellular_sodium_concentration", 155);

        // for endo, m checking
        p_model->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", 0.0225); // 0.09
        //p_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", 0.0225); // 0.09
        
        //p_model->SetParameter("ca_activated_cl_current_conductance", 1.54); 
        //p_model->SetParameter("late_sodium_current_conductance", 0.015);
/*
        p_model->SetParameter("membrane_sodium_calcium_exchanger_current_conductance", 6.2); // ncx 5.5
        p_model->SetParameter("membrane_inward_rectifier_potassium_current_conductance", 0.57153); // k1
        p_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", 0.135); // 192 // 0.152 gkr
        p_model->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", 0.022); // 0.054 gks
       // p_model->SetParameter("membrane_ca_activated_potassium_current_conductance", 0.054813);
        p_model->SetParameter("extracellular_calcium_concentration", 2.8); // 2.2  // 2.7
        p_model->SetParameter("extracellular_potassium_concentration", 5.7); // 5.9
        p_model->SetParameter("extracellular_sodium_concentration", 145); // 149
        p_model->SetParameter("membrane_inward_rectifier_potassium_current_conductance", 0.57); //0.57153

        p_model->SetParameter("membrane_slow_transient_outward_current_conductance", 0.0); // 0
        p_model->SetParameter("membrane_fast_transient_outward_current_conductance", 0.0); // 0 //0.1144
        */

       // ischemica

        bool ischemia = false, phase1a = false, phase1b = false; 

        if (ischemia == true)
        {
                double decrease_factor_gna = 0;
                double k_o = 0;
                double fkatp = 0;
                double decrease_factor_ncx = 0;
                double decrease_factor_nak = 0;
                double gnal_factor = 0;
                //double decrease_factor_iup = 0;
                //double decrease_factor_irel = 0;
                //double gcacl = 0;

                double gna = p_model->GetParameter("membrane_fast_sodium_current_conductance");
                double gca = p_model->GetParameter("membrane_L_type_calcium_current_conductance");

                if (phase1a == true)
                {
                        double gnal = p_model->GetParameter("membrane_late_sodium_current_conductance");

                        decrease_factor_gna = 0.85; //0.80; // 5po = 0.80 // 10po = 0.75 // 30po = 0.5 
                        k_o = 6.5; //6.3; //6.8; //7.6;
                        fkatp = 0.06; //0.14; //0.07;
                        //gcacl = 0.8; // 0.8 for epi, 0.2 for endo, 0.8 for m
                        //p_model->SetParameter("ca_activated_cl_current_conductance", gcacl);
                        gnal_factor = 3.5;

                        p_model->SetParameter("membrane_late_sodium_current_conductance", gnal*gnal_factor);
                }

                if (phase1b == true)
                {
                        double kncx = p_model->GetParameter("membrane_sodium_calcium_exchanger_current_conductance");
                        double pnak = p_model->GetParameter("membrane_sodium_potassium_pump_current_permeability");
                        //double srup = p_model->GetParameter("SR_uptake_current_max");
                        //double srrel = p_model->GetParameter("...");
                        //double nai = p_model->GetParameter("cytosolic_sodium_concentration");

                        decrease_factor_gna = 0.50; // 5po = 0.80 // 10po = 0.75 // 30po = 0.5   
                        k_o = 8.5; //8.7;//9.5;
                        fkatp = 0.1;
                        decrease_factor_ncx = 0.3; //0.2;
                        decrease_factor_nak = 0.3;
                        //decrease_factor_iup = 0.9;
                        //decrease_factor_irel = 0.05;

                        p_model->SetParameter("membrane_sodium_calcium_exchanger_current_conductance", kncx*decrease_factor_ncx);
                        p_model->SetParameter("membrane_sodium_potassium_pump_current_permeability", pnak*decrease_factor_nak);
                        //p_model->SetParameter("cytosolic_sodium_concentration", nai*increase_factor);
                        //p_model->SetParameter("background_ca_current_conductance", cab*increase_factor2);
                        //p_model->SetParameter("SR_calcium_release", srrel*decrease_factor5);
                        //p_model->SetParameter("SR_uptake_current_max", srup*decrease_factor4);  
                }

                p_model->SetParameter("membrane_L_type_calcium_current_conductance", gca*decrease_factor_gna);
                p_model->SetParameter("membrane_fast_sodium_current_conductance", gna*decrease_factor_gna);
                p_model->SetParameter("extracellular_potassium_concentration", k_o); // ctrl = 5.8 // 5po = 7 (7.6) // 10po = 8 (new 9) // 30po = 12 (new 13)
                p_model->SetParameter("membrane_atp_dependent_potassium_current_conductance", fkatp); // 0.07 all (new 0.1)

        }
        
        double max_timestep = 0.02; // 0.02
        double sampling_timestep = 0.02;//max_timestep;//0.1; //max_timestep; // it affects accuracy in a little rate
        //
        std::cout << "abs_t and rel_t = " << abs_t << std::endl;
        std::cout << "max_timestep: " << max_timestep << ", sampling_timestep: " << sampling_timestep << std::endl;
        std::cout << "layer: " << layer << std::endl;
        std::cout << "ischemia: " << ischemia << " phase1a: " << phase1a << " phase1b: " << phase1b << std::endl;

        SteadyStateRunner steady_runner(p_model);
        steady_runner.SetMaxNumPaces(5000u); //100000u // 5000 is optimal 
        bool result;
        result = steady_runner.RunToSteadyState();

        TS_ASSERT_EQUALS(result, true);

        p_model->SetMaxTimestep(max_timestep);

        double start_time = 0.0; // 4500 gets the same results as if start from 0
        double end_time = 400.0; // 5000

        OdeSolution solution = p_model->Compute(start_time, end_time, sampling_timestep); // Compute

        //solution.CalculateDerivedQuantitiesAndParameters(moda);

        //std::vector<double> CaL= solution.GetAnyVariable("membrane_L_type_calcium_current");
        //std::vector<double> Naf= solution.GetAnyVariable("membrane_fast_sodium_current");
        //std::vector<double> Gkr= solution.GetAnyVariable("membrane_rapid_delayed_rectifier_potassium_current");
        //std::vector<double> Gks= solution.GetAnyVariable("membrane_slow_delayed_rectifier_potassium_current");
        //std::vector<double> I_K1= solution.GetAnyVariable("membrane_inward_rectifier_potassium_current");
        //std::vector<double> fkatp= solution.GetAnyVariable("membrane_atp_dependent_potassium_current");

        //std::vector<double> CaL = solution.rGetDerivedQuantities("membrane_L_type_calcium_current");

        /*std::vector<std::string> output_variables;
        output_variables.push_back("membrane_L_type_calcium_current");
        HeartConfig::Instance()->SetOutputVariables(output_variables);
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
        HeartConfig::Instance()->SetVisualizeWithVtk(true);*/

        if (ischemia == false)
        {
                solution.WriteToFile(layer+"_ctrl", layer+"_ctrl", "ms");//, 1, false, 3, true);
        }

        if (ischemia && phase1a)
        {
                solution.WriteToFile(layer+"_phase1a", layer+"_phase1a", "ms");//, 1, false, 3, true);
        }

        if (ischemia && phase1b)
        {
                solution.WriteToFile(layer+"_phase1b", layer+"_phase1b", "ms");//, 1, false, 3, true);
        }        

        unsigned voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
        std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index);
        CellProperties cell_props(voltages, solution.rGetTimes());
        
        double apd25 = cell_props.GetLastActionPotentialDuration(25);
        double apd30 = cell_props.GetLastActionPotentialDuration(30);
        double apd50 = cell_props.GetLastActionPotentialDuration(50);
        //double apd75 = cell_props.GetLastActionPotentialDuration(75);
        double apd90 = cell_props.GetLastActionPotentialDuration(90);
        double apd95 = cell_props.GetLastActionPotentialDuration(95);
        //double apd100 = cell_props.GetLastActionPotentialDuration(100);
        double max_upstroke_velocity = cell_props.GetLastMaxUpstrokeVelocity();
        double time_max_upstroke_velocity = cell_props.GetTimeAtLastMaxUpstrokeVelocity();
        double peak_ap = cell_props.GetLastPeakPotential();
        double rest_ap = cell_props.GetLastRestingPotential();
        //double ampl_ap = cell_props.GetLastActionPotentialAmplitude();
        //std::vector<double> ccl = cell_props.GetCycleLengths(); 

        //std::cout << "cycle_length = " << ccl << " ms" << std::endl;

        std::cout << "peak_ap = " << peak_ap << " mV, tagret = < 67.5" << std::endl;
        std::cout << "rest_ap = " << rest_ap << " mV, tagret = -84.5 < -79.7" << std::endl;
        //std::cout << "ampl_ap = " << ampl_ap << " mV" << std::endl;
        std::cout << "stim_threshhold = " << peak_ap + rest_ap << " mV" << std::endl;

        std::cout << "Max upstroke velocity = " << max_upstroke_velocity << " mV/ms, , tagret = < 350" << std::endl;
        std::cout << "time at Max upstroke velocity = " << time_max_upstroke_velocity << " ms, " << std::endl;

        std::cout << "APD25 = " << apd25 << " ms, tagret = 27.4 - 36.5" << std::endl;
        std::cout << "APD30 = " << apd30 << " ms, tagret = 27.4 - 36.5" << std::endl;
        std::cout << "APD50 = " << apd50 << " ms, tagret = 177 - 186" << std::endl;
        //std::cout << "APD75 = " << apd75 << " ms, tagret = 215 - 223" << std::endl;
        std::cout << "APD90 = " << apd90 << " ms, tagret = 228 - 237" << std::endl;
        std::cout << "APD95 = " << apd95 << " ms, tagret = 234 - 244" << std::endl;
        std::cout << "APD_tri = " << apd30/apd90 << " ... " << std::endl;
        //std::cout << "APD100 = " << apd100 << " ms, tagret = 245 - 255" << std::endl;
        
        
        //std::cout << "K_o = " << apd25 << " ms, tagret = 27.4 - 36.5" << std::endl;
        

        //TS_ASSERT_DELTA(apd, 212.411, 1e-2);
        //TS_ASSERT_DELTA(upstroke_velocity, 338.704, 1.25);

#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};