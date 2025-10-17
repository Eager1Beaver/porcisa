#ifndef TEST_CONVERGENCE_HPP_
#define TEST_CONVERGENCE_HPP_

/*
    Mesh-resolution convergence test on a small ischemic tissue block.
    
    What it does
    - Builds a regular 3D slab mesh at several spatial steps (mm)
     - Launches a planar-ish wave using a small spherical stimulus region
    - Records activation time (AT; first V-crossing at -30 mV) at probe points
    - Uses the finest grid as the reference to compute per-probe AT error
    - Writes tidy CSV: per-resolution ATs + summary (L2 / Linf errors)

    Notes
    - Defaults to Monodomain for speed
 */

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "CommandLineArguments.hpp"

#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "AbstractCvodeCell.hpp"
#include "SimpleStimulus.hpp"
#include "UblasIncludes.hpp"

#include "DistributedTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp"

#include "pig_ventr_ap_endoCvodeOpt.hpp"

// ----------------------------- Config -----------------------------
struct ConvergenceConfig
{
    std::string outdir = "output/convergence";
    // mm steps. The last (smallest) is the reference.
    std::vector<double> res_mm = {1.0, 0.8, 0.6, 0.4, 0.2};

    // Block geometry (mm)
    double size_x_mm = 10.0; 
    double size_y_mm = 10.0;
    double size_z_mm = 10.0;

    // Time (ms)
    double t_end_ms  = 25.0;
    double dt_ode_ms = 0.02;
    double dt_pde_ms = 0.1;
    double dt_out_ms = 0.2;

    // Stimulus: spherical region in one corner to create a wavefront
    double stim_amp_uAcm3 = -45000.0;
    double stim_dur_ms    = 2.0;
    c_vector<double,3> stim_center = {{0.8, 0.8, 0.8}}; // mm from origin
    double stim_radius_mm = 0.6;

    // Activation-time threshold (mV) for first upstroke crossing
    double at_threshold_mV = -30.0;

    // Probe locations (mm) — by default 5 points along the main diagonal
    std::vector<c_vector<double,3>> probes_mm = {
        {{2.0, 2.0, 2.0}},
        {{4.0, 4.0, 4.0}},
        {{6.0, 6.0, 6.0}},
        {{8.0, 8.0, 8.0}},
        {{9.0, 9.0, 9.0}}
    };
};

// ----------------------------- Cell factory -----------------------------
class HomogPorcineFactory : public AbstractCardiacCellFactory<3>
{
public:
    HomogPorcineFactory(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                        boost::shared_ptr<AbstractStimulusFunction> pStim)
    : AbstractCardiacCellFactory<3>(pSolver), mpStim(pStim) {}

    AbstractCardiacCellInterface* CreateCardiacCellForTissueNode(Node<3>* pNode) override
    {
        // All cells identical (endo) for a clean, homogeneous AT field
        return new pig_ventr_ap_endoCvodeOpt(this->mpSolver, mpStim);
    }

private:
    boost::shared_ptr<AbstractStimulusFunction> mpStim;
};

// ----------------------------- Helpers -----------------------------
inline std::string ToCsv(const std::vector<double>& v)
{
    std::ostringstream ss;
    for (size_t i=0; i<v.size(); ++i) { if (i) ss << ","; ss << std::setprecision(10) << v[i]; }
    return ss.str();
}

inline double L2Error(const std::vector<double>& ref, const std::vector<double>& x)
{
    double s=0; for (size_t i=0;i<ref.size();++i){ double d=x[i]-ref[i]; s+=d*d; }
    return std::sqrt(s/static_cast<double>(ref.size()));
}
inline double LinfError(const std::vector<double>& ref, const std::vector<double>& x)
{
    double m=0; for (size_t i=0;i<ref.size();++i){ m=std::max(m, std::fabs(x[i]-ref[i])); }
    return m;
}

inline unsigned NearestNodeIndex(const DistributedTetrahedralMesh<3,3>& mesh,
                                 const c_vector<double,3>& xyz)
{
    // Simple linear search
    double best = 1e300; unsigned arg=0;
    for (unsigned i=0; i<mesh.GetNumNodes(); ++i)
    {
        const c_vector<double,3>& p = mesh.GetNode(i)->rGetLocation();
        const double dx=p[0]-xyz[0], dy=p[1]-xyz[1], dz=p[2]-xyz[2];
        const double d2=dx*dx+dy*dy+dz*dz;
        if (d2<best){ best=d2; arg=i; }
    }
    return arg;
}

inline void WriteCsv(const std::string& path, const std::string& text)
{
    std::ofstream f(path.c_str(), std::ios::out|std::ios::trunc);
    if(!f.is_open()) throw std::runtime_error("Cannot open CSV for writing: "+path);
    f << text;
}

// ----------------------------- Core runner -----------------------------
/*
    Run the convergence test and write CSVs:
        - AT_per_resolution.csv : header "res_mm,AT_probe0,...,AT_probeN"
        - AT_error_summary.csv  : header "res_mm,L2,Linf"
 */
inline void RunConvergence(const ConvergenceConfig& cfg)
{
    // Make sure resolutions are sorted and unique
    auto res = cfg.res_mm;
    std::sort(res.begin(), res.end());
    res.erase(std::unique(res.begin(), res.end()), res.end());
    if (res.size() < 2) throw std::runtime_error("Need at least two resolutions for convergence.");

    // Prepare output
    OutputFileHandler handler(cfg.outdir, /*cleanOutputDirectory=*/true);
    const std::string path_at  = handler.GetOutputDirectoryFullPath() + "AT_per_resolution.csv";
    const std::string path_err = handler.GetOutputDirectoryFullPath() + "AT_error_summary.csv";

    // Storage of ATs per resolution
    std::vector<std::vector<double>> ats_all;
    ats_all.reserve(res.size());

    // Main loop over resolutions
    for (double h_mm : res)
    {
        if (PetscTools::AmMaster())
            std::cout << "[Convergence] Resolving at h=" << h_mm << " mm\n";

        // Build a regular slab mesh [0,Sx]x[0,Sy]x[0,Sz] with step h
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructRegularSlabMesh(h_mm / 1000.0,     // spacing in meters
                                      cfg.size_x_mm/1000.0,
                                      cfg.size_y_mm/1000.0,
                                      cfg.size_z_mm/1000.0);

        // Build a spherical stimulus

        // Precompute which nodes are stimulated
        std::vector<char> stim_mask(mesh.GetNumNodes(), 0);
        const double r2 = cfg.stim_radius_mm*cfg.stim_radius_mm;
        for (unsigned i=0; i<mesh.GetNumNodes(); ++i)
        {
            const c_vector<double,3>& p = mesh.GetNode(i)->rGetLocation();
            c_vector<double,3> pm = p; // meters -> convert to mm for comparison
            pm[0]*=1000.0; pm[1]*=1000.0; pm[2]*=1000.0;
            const double dx=pm[0]-cfg.stim_center[0], dy=pm[1]-cfg.stim_center[1], dz=pm[2]-cfg.stim_center[2];
            if (dx*dx+dy*dy+dz*dz <= r2) stim_mask[i]=1;
        }

        // Build per-node stimulus functions
        std::vector<boost::shared_ptr<SimpleStimulus>> node_stims(mesh.GetNumNodes());
        for (unsigned i=0; i<mesh.GetNumNodes(); ++i)
        {
            double amp = stim_mask[i] ? cfg.stim_amp_uAcm3 : 0.0;
            node_stims[i].reset(new SimpleStimulus(amp, cfg.stim_dur_ms, 0.5 /*start*/));
        }

        // Factory that uses the per-node stimulus
        class NodeStimFactory : public AbstractCardiacCellFactory<3>
        {
        public:
            NodeStimFactory(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                            const std::vector<boost::shared_ptr<SimpleStimulus>>* pStims)
            : AbstractCardiacCellFactory<3>(pSolver), mpStims(pStims) {}
            AbstractCardiacCellInterface* CreateCardiacCellForTissueNode(Node<3>* pNode) override
            {
                return new pig_ventr_ap_endoCvodeOpt(this->mpSolver, (*mpStims)[pNode->GetIndex()]);
            }
        private:
            const std::vector<boost::shared_ptr<SimpleStimulus>>* mpStims;
        } factory(boost::shared_ptr<AbstractIvpOdeSolver>(), &node_stims);

        // Problem
        MonodomainProblem<3> problem(&factory);
        problem.SetMesh(&mesh);

        // Timesteps
        HeartConfig::Instance()->SetOdeTimeStep(cfg.dt_ode_ms);
        HeartConfig::Instance()->SetPdeTimeStep(cfg.dt_pde_ms);
        HeartConfig::Instance()->SetPrintingTimeStep(cfg.dt_out_ms);
        HeartConfig::Instance()->SetSimulationDuration(cfg.t_end_ms);

        // Initialise and solve
        problem.Initialise();
        problem.Solve();

        // Build probe node index list for this mesh
        std::vector<unsigned> probe_nodes(cfg.probes_mm.size());
        for (size_t k=0; k<cfg.probes_mm.size(); ++k)
        {
            c_vector<double,3> p = cfg.probes_mm[k];
            // convert mm → m for mesh coords
            p[0]/=1000.0; p[1]/=1000.0; p[2]/=1000.0;
            probe_nodes[k] = NearestNodeIndex(mesh, p);
        }

        // Extract ATs from voltage files: Chaste writes Vm nodes over time; read from memory.
        // Re-simulate a short sampling at dt_out and 
        // track first threshold crossing for the probe nodes using the cell states.

        std::vector<double> at_ms(cfg.probes_mm.size(), std::numeric_limits<double>::quiet_NaN());
        const double thr = cfg.at_threshold_mV;

        // Time loop
        {
            // Reset problem to t=0 and step again for sampling;
            MonodomainProblem<3> sampler(&factory);
            sampler.SetMesh(&mesh);
            HeartConfig::Instance()->SetSimulationDuration(cfg.t_end_ms);
            HeartConfig::Instance()->SetOdeTimeStep(cfg.dt_ode_ms);
            HeartConfig::Instance()->SetPdeTimeStep(cfg.dt_pde_ms);
            HeartConfig::Instance()->SetPrintingTimeStep(cfg.dt_out_ms);
            sampler.Initialise();

            double t = 0.0;
            const unsigned nsteps = static_cast<unsigned>(std::ceil(cfg.t_end_ms / cfg.dt_out_ms));
            for (unsigned s=0; s<=nsteps; ++s)
            {
                // Access Vm for probes and detect first upward crossing
                for (size_t k=0; k<probe_nodes.size(); ++k)
                {
                    if (!std::isnan(at_ms[k]))
                        continue; // already found

                    AbstractCardiacCellInterface* p_cell =
                        sampler.rGetTissue()->GetCardiacCell(probe_nodes[k]);
                    const double vm = p_cell->GetVoltage();
                    // We need previous sample to detect crossing; 
                    // For simplicity, treat <thr at previous step and >=thr now as crossing.
                    // Initialize prev below threshold.
                    static std::vector<double> prev;
                    if (prev.size()!=probe_nodes.size()) prev.assign(probe_nodes.size(), thr-100.0);
                    if (prev[k] < thr && vm >= thr) at_ms[k] = t;
                    prev[k] = vm;
                }

                if (std::all_of(at_ms.begin(), at_ms.end(), [](double v){return !std::isnan(v);} ))
                    break;

                sampler.SolveOneTimeStep();
                t += cfg.dt_out_ms;
            }
        }

        ats_all.push_back(at_ms);
    }

    // Reference = last (finest) set
    const std::vector<double>& at_ref = ats_all.back();

    // Write ATs per resolution
    {
        std::ostringstream csv;
        csv << "res_mm";
        for (size_t k=0; k<at_ref.size(); ++k) csv << ",AT_probe" << k;
        csv << "\n";
        for (size_t i=0; i<res.size(); ++i)
        {
            csv << std::fixed << std::setprecision(3) << res[i] << ",";
            csv << ToCsv(ats_all[i]) << "\n";
        }
        if (PetscTools::AmMaster()) WriteCsv(path_at, csv.str());
    }

    // Error summary vs reference
    {
        std::ostringstream csv;
        csv << "res_mm,L2,Linf\n";
        for (size_t i=0; i<res.size(); ++i)
        {
            const double l2   = L2Error(at_ref, ats_all[i]);
            const double linf = LinfError(at_ref, ats_all[i]);
            csv << std::fixed << std::setprecision(3) << res[i] << ","
                << std::setprecision(10) << l2 << "," << linf << "\n";
        }
        if (PetscTools::AmMaster()) WriteCsv(path_err, csv.str());
    }

    if (PetscTools::AmMaster())
        std::cout << "[Convergence] Wrote:\n  - " << path_at << "\n  - " << path_err << std::endl;
}

// ----------------------------- Main -----------------------------
#ifdef PORCISA_CONVERGENCE_MAIN
int main(int argc, char* argv[])
{
    CommandLineArguments* p_args = CommandLineArguments::Instance();
    p_args->Parse(argc, argv);

    ConvergenceConfig cfg;
    if (p_args->OptionExists("--out"))  cfg.outdir = p_args->GetOption("--out");
    if (p_args->OptionExists("--tend")) cfg.t_end_ms = atof(p_args->GetOption("--tend").c_str());
    if (p_args->OptionExists("--res"))
    {
        // --res 1.0,0.8,0.6,0.4,0.2
        std::string s = p_args->GetOption("--res");
        std::replace(s.begin(), s.end(), ';', ',');
        std::stringstream ss(s);
        cfg.res_mm.clear();
        for (std::string tok; std::getline(ss, tok, ','); )
            if (!tok.empty()) cfg.res_mm.push_back(atof(tok.c_str()));
    }

    try {
        RunConvergence(cfg);
    } catch (const std::exception& e) {
        if (PetscTools::AmMaster())
            std::cerr << "Fatal: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
#endif // PORCISA_CONVERGENCE_MAIN

#endif // TEST_CONVERGENCE_HPP_
