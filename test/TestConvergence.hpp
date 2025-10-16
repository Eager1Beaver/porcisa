#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "TetrahedralMesh.hpp"
#include "SimpleStimulus.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "carro_2011_epi_pig_scratch_4_go1CvodeOpt.hpp"
//#include "carro_2011_endo_pig_scratch_4_go1CvodeOpt.hpp"
//#include "carro_2011_m_pig_scratch_4_go1CvodeOpt.hpp"

#ifdef CHASTE_CVODE

class BenchmarkCellFactory : public AbstractCardiacCellFactory<3> // <3> here
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

    bool mIschaemia;
    
    double  _k0;
    double  _ikatp;
    double  _ina;
    double  _ica;
    double  _gncx;
    double  _pnak;
    double  _inal;

public:
    BenchmarkCellFactory(   bool ischaemia,
                            double  k0,
                            double  ikatp,
                            double  ina,
                            double  ica,
                            double  gncx,
                            double  pnak,
                            double  inal
                            )
        : AbstractCardiacCellFactory<3>(), // <3> here as well!
          mpStimulus(new SimpleStimulus(-100000.0, 2)),
          _k0(k0),
        _ikatp(ikatp),
        _ina(ina),
        _ica(ica),
        _gncx(gncx),
        _pnak(pnak),
        _inal(inal)
    {
        mIschaemia = ischaemia;
    }

    AbstractCvodeCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        AbstractCvodeCell* p_cell;
        boost::shared_ptr<AbstractIvpOdeSolver> p_empty_solver;

        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        double z = pNode->rGetLocation()[2];

        if ((x<0.1+1e-6) && (y<0.1+1e-6) && (z<0.1+1e-6))
        {
            p_cell = new Cellcarro_2011_epi_pig_scratch_4_go1FromCellMLCvodeOpt(p_empty_solver, mpStimulus);
        }
        else
        {
            p_cell = new Cellcarro_2011_epi_pig_scratch_4_go1FromCellMLCvodeOpt(p_empty_solver, mpZeroStimulus);
        }

        double k_i = _k0;
        double fkatp_i = _ikatp;
        double ina_i = _ina;
        double ica_i = _ica;
        double gncx_i = _gncx;
        double pnak_i = _pnak;
        double inal_i = _inal;

        if(mIschaemia)
                {
                    // ICZ
                    p_cell->SetParameter("extracellular_potassium_concentration",k_i);
                    double gna = p_cell->GetParameter("membrane_fast_sodium_current_conductance");
                    double gca = p_cell->GetParameter("membrane_L_type_calcium_current_conductance");
                    p_cell->SetParameter("membrane_fast_sodium_current_conductance",gna*ina_i);
                    p_cell->SetParameter("membrane_L_type_calcium_current_conductance",gca*ica_i);
                    p_cell->SetParameter("membrane_atp_dependent_potassium_current_conductance",fkatp_i);

                    double get_gNCX = p_cell->GetParameter("membrane_sodium_calcium_exchanger_current_conductance");
                    p_cell->SetParameter("membrane_sodium_calcium_exchanger_current_conductance", get_gNCX * gncx_i);
                    double get_pNaK = p_cell->GetParameter("membrane_sodium_potassium_pump_current_permeability");
                    p_cell->SetParameter("membrane_sodium_potassium_pump_current_permeability", get_pNaK * pnak_i);

                    double get_gNaL = p_cell->GetParameter("membrane_late_sodium_current_conductance");
                    p_cell->SetParameter("membrane_late_sodium_current_conductance", get_gNaL * inal_i);
                 }

        p_cell->SetTolerances(1e-5,1e-7);

        return p_cell;
    }
};

#endif // CHASTE_CVODE

class TestConvergence : public CxxTest::TestSuite
{
public:
    void TestConvergence()
    {
#ifdef CHASTE_CVODE

        bool ischaemia = true;

        double k0       = 0.0;
        double ikatp    = 0.0;
        double ina      = 0.0;
        double ica      = 0.0;
        double gncx     = 0.0;    
        double pnak     = 0.0;     
        double inal     = 0.0;     
        double factor_h = 1.0; 

        if (ischaemia == true) 
            {
                k0 =    8.3; // default 8.3
                ikatp = 0.1; // default 0.1 
                ina =   0.5; // default 0.5 
                ica =   0.5; // default 0.5 
                inal =  3.0; // default 3.0
                gncx =  0.33; // default 0.33 
                pnak =  0.33; // default 0.33

                factor_h = 0.1;
            }

            double intra_longi = 3.35; double intra_trans = 1.755; double intra_normal = 0.8775;
            double extra_longi = 12.4; double extra_trans = 5.5; double extra_normal = 5.5;
            HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(intra_longi*factor_h, intra_trans*factor_h, intra_normal*factor_h));
            HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(extra_longi, extra_trans, extra_normal));


        DistributedTetrahedralMesh<3,3> mesh;
        double h=0.02; // 0.02 cm = 200 um
        mesh.ConstructRegularSlabMesh(h, 1 /*length*/, 1 /*width*/, 1 /*depth*/);
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);

        HeartConfig::Instance()->SetSimulationDuration(30); //ms
        HeartConfig::Instance()->SetOutputDirectory("TestConvIsch_200um");
        HeartConfig::Instance()->SetOutputFilenamePrefix("conv");

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.02, 0.02, 1);

        HeartConfig::Instance()->SetVisualizeWithVtk(true);
	    HeartConfig::Instance()->SetVisualizeWithParallelVtk(false);

        // Activation time map
        std::vector<double> upstroke_time_map;
        upstroke_time_map.push_back(0.0);                          
        HeartConfig::Instance()->SetUpstrokeTimeMaps(upstroke_time_map);
        
        // Max Upstroke Velocity map
        std::vector<double> maxupstroke_vel_map;
        maxupstroke_vel_map.push_back(0.0);                             
        HeartConfig::Instance()->SetMaxUpstrokeVelocityMaps(maxupstroke_vel_map);
        
        // CV map
        std::vector<unsigned> cv_map;
        double  ax= 0.01 ,   ay= 0.01,   az= 0.01;
        ChastePoint<3> cv_point(ax,ay,az);
        int cv_point_idx;
        cv_point_idx = mesh.GetNearestNodeIndex(cv_point);
        cv_map.push_back(cv_point_idx);                     
        HeartConfig::Instance()->SetConductionVelocityMaps(cv_map);

        std::string preconditioner;
        std::string solver;    
        //preconditioner = "blockdiagonal";
        solver = "symmlq";
        //HeartConfig::Instance()->SetKSPPreconditioner( preconditioner.c_str() );
        HeartConfig::Instance()->SetKSPSolver( solver.c_str() );

        BenchmarkCellFactory cell_factory(  ischaemia,
                                            k0,
                                            ikatp,
                                            ina,
                                            ica,
                                            gncx,
                                            pnak,
                                            inal
                                            );

        MonodomainProblem<3> monodomain_problem( &cell_factory );
        monodomain_problem.SetMesh( &mesh );

        std::cout << "ischaemia: " << ischaemia << std::endl;

        monodomain_problem.SetWriteInfo();
        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();

#else
        std::cout << "CVODE is not installed, or CHASTE is not configured to use it, check your hostconfig settings." << std::endl;
        // TS_ASSERT(false); // uncomment if you want to ensure CVODE is set up on your system.
#endif // CHASTE_CVODE
    }
};