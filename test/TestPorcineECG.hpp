#ifndef TEST_PORCINE_ECG_HPP_
#define TEST_PORCINE_ECG_HPP_

/**
 	Run porcine heart-torso simulations with 
		heterogeneity, ischemia geometry, and multi-focus root stimuli.
 
 	Control at the top (ExperimentConfig):
   - Output directory & run time
   - Mesh filename (prefix; no extension)
   - Heterogeneity toggles + apex/base anchors
   - Ischemia: core/border radii and center
   - Stimuli: spherical roots (radius 0.25) & default start times
   - Optional: ECG electrodes (surface leads)
 */

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <fstream>
#include <cctype>

#include "HeartConfig.hpp"
#include "BidomainProblem.hpp"
#include "CommandLineArguments.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"

#include "PorcineHeterogen.hpp"
#include "PorcineConductivityModifier.hpp"
#include "PorcineCellFactory.hpp"

// ------------ configuration ------------
struct ExperimentConfig
{
    // I/O
    std::string output_dir      = "output/Run_01";
    std::string mesh_prefix     = "mesh/pig_mesh_400k";

    // Timing (ms)
    double simulation_duration_ms = 300.0;  // total time
    double printing_dt_ms         = 1.0;    // sampling/printing
    double ode_dt_ms              = 0.1;   // ODE timestep
    double pde_dt_ms              = 0.1;    // PDE timestep

    // Heterogeneity
    bool enable_apex_base         = true;
    bool enable_transmural        = false;  // stays neutral unless you wire endo/epi distances
    bool enable_interventricular  = true;
    double apex_z                 = 47.246;   // apex coordinate in mesh space
    double base_z                 = 54.383;   // base coordinate in mesh space

    // Ischemia geometry (sphere in heart coords)
    bool   ischaemia_on           = true;
    double ischemia_cx            = 23.0;
    double ischemia_cy            = 4.9;
    double ischemia_cz            = 49.5;
    double ischemia_core_radius   = 2.3;    // ICZ
    double ischemia_border_radius = 2.8;    // BZ (>= core)

    // Ischemia domain gate 0 -> (0.1,0.5), 1 -> (>0.1), 2 -> (0.1,0.3)
    unsigned transmural_gate      = 0;

    // Stimuli
    bool use_user_root_stimulus_preset = true; // 8 roots, radius 0.25, representative starts
    double stimulus_amplitude_uAcm3     = -45000.0;
    double stimulus_duration_ms         = 2.0;

    // ECG electrodes (optional; set empty to skip)
    std::vector<c_vector<double,3>> electrodes; // torso positions, in mesh space

	// ECG via node index file (one index per line)
    // If empty, ECG electrodes are not set from a file.
    std::string ecg_node_index_file = ""; // e.g., "mesh/ecg_nodes.txt"


    // Conductivity modifier
    double conductivity_scale_in_core   = 0.6; // multiply tensor inside ICZ (e.g., 0.6)
};

// ------------ helpers ------------
inline std::string TimestampTag()
{
    std::stringstream ss;
    ss << "t" << std::setfill('0') << std::setw(2) << PetscTools::GetMyRank();
    return ss.str();
}

inline void ConfigureHeart(const ExperimentConfig& cfg)
{
    HeartConfig* hc = HeartConfig::Instance();

    // Mesh & problem type
    hc->SetMeshFileName(cfg.mesh_prefix, cp::media_type::Orthotropic); // cp::media_type::Orthotropic)
    hc->SetSimulationDuration(cfg.simulation_duration_ms);
    hc->SetOutputDirectory(cfg.output_dir);

    // Timesteps
    hc->SetPdeTimeStep(cfg.pde_dt_ms);
    hc->SetOdeTimeStep(cfg.ode_dt_ms);
    hc->SetPrintingTimeStep(cfg.printing_dt_ms);

    // Conductivities
    // hc->SetConductivity( tissue, values... ) // (optional)

    // Surface ECG (optional)
    if (!cfg.electrodes.empty())
    {
        for (auto& e : cfg.electrodes)
        {
            hc->SetElectrode(e);
        }
        hc->SetUseSurface(ecg::yes);
    }
    else
    {
        // If computing pseudo-ECG between two points:
        // hc->SetElectrode(...); hc->SetUsePseudoEcg(true);
    }

    // Solver & writing options
    hc->SetVisualizeWithMeshalyzer(false);
    hc->SetVisualizeWithVtk(true);
    hc->SetWriteEcgPdfs(false);
    hc->SetCheckCompatibility(false);
}

//
inline std::vector<unsigned> LoadNodeIndicesFile(const std::string& path)
{
    std::vector<unsigned> idxs;
    if (path.empty()) return idxs;

    std::ifstream fin(path.c_str());
    if (!fin.is_open())
        throw std::runtime_error("Cannot open ECG node index file: " + path);

    std::string line;
    while (std::getline(fin, line))
    {
        // strip comments
        auto hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);

        // trim
        auto notspace = [](int ch){ return !std::isspace(ch); };
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), notspace));
        line.erase(std::find_if(line.rbegin(), line.rend(), notspace).base(), line.end());

        if (line.empty()) continue;

        // parse unsigned
        unsigned v = 0;
        std::stringstream ss(line);
        ss >> v;
        if (!ss.fail()) idxs.push_back(v);
    }
    return idxs;
}


// ------------ main ------------
inline void RunPorcineExperiment(const ExperimentConfig& cfg)
{
    ConfigureHeart(cfg);

    // Build problem
    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
    HeartGeometryInformation<3>* p_geom = HeartConfig::Instance()->GetHeartGeometryInformation();

    // Heterogeneity
    PorcineHeterogen hetero(p_geom,
                            cfg.enable_apex_base,
                            cfg.enable_transmural,
                            cfg.enable_interventricular);
    if (cfg.enable_apex_base)
    {
        hetero.SetApexBaseAnchors(cfg.apex_z, cfg.base_z);
    }

    // Cell factory
    PorcineCellFactory factory(p_solver, &hetero, cfg.ischaemia_on);
    factory.SetStimulusAmplitude(cfg.stimulus_amplitude_uAcm3);
    factory.SetStimulusDurationMs(cfg.stimulus_duration_ms);
    if (cfg.use_user_root_stimulus_preset)
    {
        factory.UseUserRootStimulusPreset();
    }

    // Ischemia geometry & transmural gating (core/border)
    factory.SetIschemiaGeometry(cfg.ischemia_cx, cfg.ischemia_cy, cfg.ischemia_cz,
                                cfg.ischemia_core_radius, cfg.ischemia_border_radius);
    factory.SetTransmuralMode(cfg.transmural_gate);

    // Optional: adjust ischemia multipliers
    // factory.SetIschemiaMultipliersCore( 1.2, 0.8, 0.8, 1.5, 1.2, 0.7, 1.5 );
    // factory.SetIschemiaMultipliersBorder(1.1, 0.9, 0.9, 1.3, 1.1, 0.9, 1.2 );

    // Problem
    BidomainProblem<3> problem(&factory);

    // Initialisation builds the tissue; after this, you can attach conductivity modifiers that
    // need access to the fully constructed problem/mesh.
    problem.Initialise();

	// Wire geometry into heterogeneity
	hetero.SetGeometry(factory.GetHeartGeometryInformation());


	    // ----- ECG electrodes from node index file (optional) -----
    if (!cfg.ecg_node_index_file.empty())
    {
        const std::vector<unsigned> ecg_nodes = LoadNodeIndicesFile(cfg.ecg_node_index_file);
        if (ecg_nodes.empty())
        {
            if (PetscTools::AmMaster())
                std::cout << "ECG file provided but no indices parsed: "
                          << cfg.ecg_node_index_file << std::endl;
        }
        else
        {
            // Alternative options
            // 1) problem.rGetEcgCalculator()->AddElectrodeAtNode(node_index);
            // 2) problem.AddSurfaceElectrodeAtNode(node_index);
            // 3) problem.rGetTissue()->AddSurfaceElectrodeAtNode(node_index);
            //
            bool attached_by_node = false;
            // for (unsigned ni : ecg_nodes)
            // {
            //     problem.rGetEcgCalculator()->AddElectrodeAtNode(ni);
            //     attached_by_node = true;
            // }

            if (!attached_by_node)
            {
                // Convert node indices to coordinates and attach by position.
                for (unsigned ni : ecg_nodes)
                {
                    const c_vector<double,3>& xyz = problem.rGetMesh().GetNode(ni)->rGetLocation();
                    // Options:
                    // a) problem.rGetEcgCalculator()->AddSurfaceElectrode(xyz);
                    // b) HeartConfig::Instance()->SetElectrode(xyz);
                    //
                    // (a) goes if post-initialisation:
                    // problem.rGetEcgCalculator()->AddSurfaceElectrode(xyz);
                    HeartConfig::Instance()->SetElectrode(xyz);
                }
                HeartConfig::Instance()->SetUseSurface(ecg::yes);
            }

            if (PetscTools::AmMaster())
            {
                std::cout << "Registered " << ecg_nodes.size()
                          << " ECG electrodes from " << cfg.ecg_node_index_file << "\n";
            }
        }
    }

    // Conductivity modifier (ischemic sphere)
    {
        PorcineConductivityModifier cond_mod(
            /*ischaemia*/ cfg.ischaemia_on,
            /*factor_h */ cfg.conductivity_scale_in_core,
            /*radius   */ cfg.ischemia_core_radius,
            /*cx,cy,cz */ cfg.ischemia_cx, cfg.ischemia_cy, cfg.ischemia_cz,
            /*gate     */ cfg.transmural_gate,
            /*problem  */ problem);

        // Options:
        // 1)* problem.rGetTissue()->SetConductivityModifier(&cond_mod);
        // 2) problem.SetConductivityModifier(&cond_mod);
        // 3) HeartConfig::Instance()->SetConductivityModifier(&cond_mod); // rarer

        // Solve with modifier in place
        problem.Solve();
    }

    if (PetscTools::AmMaster())
    {
        std::cout << "Finished simulation. Output: " << cfg.output_dir << std::endl;
    }
}

// ------------ optional CLI wrapper ------------
/*
	If standalone app, define PORCISA_STANDALONE_MAIN and compile this header.
 	Then run like, e.g.,
		--out output/Run02 --tend 150 --ab 52.2,55.8 --isz 23.0,6.9,52.5 --icz 1.8 --bz 2.8 --ischemia 1
 */
#ifdef PORCISA_STANDALONE_MAIN
int main(int argc, char* argv[])
{
    CommandLineArguments* p_args = CommandLineArguments::Instance();
    p_args->Parse(argc, argv);

    ExperimentConfig cfg;

    // Flag parsing
    if (p_args->OptionExists("--out"))   { cfg.output_dir = p_args->GetOption("--out"); }
    if (p_args->OptionExists("--mesh"))  { cfg.mesh_prefix = p_args->GetOption("--mesh"); }
    if (p_args->OptionExists("--tend"))  { cfg.simulation_duration_ms = atof(p_args->GetOption("--tend").c_str()); }
    if (p_args->OptionExists("--dt"))    { cfg.printing_dt_ms = atof(p_args->GetOption("--dt").c_str()); }
    if (p_args->OptionExists("--ode"))   { cfg.ode_dt_ms = atof(p_args->GetOption("--ode").c_str()); }
    if (p_args->OptionExists("--pde"))   { cfg.pde_dt_ms = atof(p_args->GetOption("--pde").c_str()); }
    if (p_args->OptionExists("--ischemia")) { cfg.ischaemia_on = atoi(p_args->GetOption("--ischemia").c_str()) != 0; }
    if (p_args->OptionExists("--icz"))   { cfg.ischemia_core_radius = atof(p_args->GetOption("--icz").c_str()); }
    if (p_args->OptionExists("--bz"))    { cfg.ischemia_border_radius = atof(p_args->GetOption("--bz").c_str()); }
    if (p_args->OptionExists("--isz"))
    {
        // cx,cy,cz
        std::string s = p_args->GetOption("--isz");
        double cx,cy,cz; char c1,c2;
        std::stringstream ss(s); ss >> cx >> c1 >> cy >> c2 >> cz;
        if (ss) { cfg.ischemia_cx = cx; cfg.ischemia_cy = cy; cfg.ischemia_cz = cz; }
    }
    if (p_args->OptionExists("--ab"))
    {
        // apex_z,base_z
        std::string s = p_args->GetOption("--ab");
        double az,bz; char c;
        std::stringstream ss(s); ss >> az >> c >> bz;
        if (ss) { cfg.apex_z = az; cfg.base_z = bz; }
    }
    if (p_args->OptionExists("--gate"))  { cfg.transmural_gate = (unsigned)atoi(p_args->GetOption("--gate").c_str()); }
    if (p_args->OptionExists("--gcond")) { cfg.conductivity_scale_in_core = atof(p_args->GetOption("--gcond").c_str()); }

	if (p_args->OptionExists("--ecg-nodes"))
    {
        cfg.ecg_node_index_file = p_args->GetOption("--ecg-nodes");
    }

	if (PetscTools::AmMaster()) {
    	std::cout << "[Porcisa] Out=" << cfg.output_dir
              << "  Tend=" << cfg.simulation_duration_ms << " ms"
              << "  Mesh=" << cfg.mesh_prefix
              << "  Ischemia=" << (cfg.ischaemia_on?"on":"off")
              << "  ICZ=" << cfg.ischemia_core_radius
              << "  BZ="  << cfg.ischemia_border_radius
              << std::endl;
			}


    try {
        RunPorcineExperiment(cfg);
    } catch (std::exception& e) {
        std::cerr << "Fatal: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
#endif // PORCISA_STANDALONE_MAIN

#endif // TEST_PORCINE_ECG_HPP_
