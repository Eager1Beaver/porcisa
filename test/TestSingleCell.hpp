#ifndef TEST_SINGLE_CELL_HPP_
#define TEST_SINGLE_CELL_HPP_

/*
        Pace a single porcine AP model to steady state and compute biomarkers.

        Outputs (CSV in cfg.outdir):
        - singlecell_biomarkers.csv   // one row with computed metrics
        - singlecell_trace.csv        // time, Vm for the final (analyzed) beat
 
        How it works
        1) Paces the cell with a RegularStimulus for up to max_beats.
        2) After each beat (period), computes APD90; checks convergence on the last K beats.
        3) Uses the final beat to compute biomarkers (RMP, Vpeak, APA, dVdt_max, APD50, APD90, triangulation).
 */

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "OutputFileHandler.hpp"
#include "CommandLineArguments.hpp"
#include "PetscTools.hpp"

#include "AbstractCvodeCell.hpp"
#include "AbstractCardiacCellInterface.hpp"
#include "RegularStimulus.hpp"
#include "SimpleStimulus.hpp"
#include "UblasIncludes.hpp"

// If switching a layer: include epi/m variants and swap the constructor.
#include "pig_ventr_ap_endoCvodeOpt.hpp"

struct SingleCellConfig
{
    std::string outdir = "Porcisa/SingleCell";
    // Pacing
    double period_ms        = 1000.0;   // BCL
    double stim_start_ms    = 50.0;     // first pulse start
    double stim_dur_ms      = 2.0;      // ms
    double stim_amp_uAcm2   = -52.0;    // uA/cm^2 (single-cell units)
    // Simulation/sampling
    unsigned max_beats      = 300;
    unsigned min_beats      = 50;
    unsigned steady_window  = 5;        // last K beats to check convergence
    double   apd_tol_ms     = 0.2;      // APD90 change tolerance for convergence
    double   dt_sample_ms   = 0.02;     // sampling step for recorded Vm
    // Upstroke detection/APD thresholds
    double   dvdt_threshold = 10.0;     // mV/ms for upstroke detection
    // Biomarker trace export
    bool     write_trace_csv = true;
};

// --------- helpers ----------
inline void WriteCsv(const std::string& path, const std::string& text)
{
    std::ofstream f(path.c_str(), std::ios::out | std::ios::trunc);
    if(!f.is_open()) throw std::runtime_error("Cannot open for writing: "+path);
    f << text;
}

struct BeatMetrics
{
    double rmp_mV = std::numeric_limits<double>::quiet_NaN();
    double vpeak_mV = std::numeric_limits<double>::quiet_NaN();
    double apa_mV = std::numeric_limits<double>::quiet_NaN();
    double dvdt_max = std::numeric_limits<double>::quiet_NaN(); // mV/ms
    double apd50_ms = std::numeric_limits<double>::quiet_NaN();
    double apd90_ms = std::numeric_limits<double>::quiet_NaN();
    double triangulation_ms = std::numeric_limits<double>::quiet_NaN(); // APD90-APD50
    double t_up_ms = std::numeric_limits<double>::quiet_NaN(); // upstroke time
};

/*
        Extract biomarkers from a single-beat trace.
        t, v must cover exactly one beat window [t0, t0+period].
 */
inline BeatMetrics AnalyzeBeat(const std::vector<double>& t,
                               const std::vector<double>& v,
                               double dvdt_threshold /*mV/ms*/)
{
    BeatMetrics bm;
    const size_t n = t.size();
    if (n < 3) return bm;

    // Finite-difference dV/dt
    std::vector<double> dvdt(n, 0.0);
    for (size_t i=1; i+1<n; ++i)
    {
        const double dt = std::max(1e-12, t[i+1]-t[i-1]);
        dvdt[i] = (v[i+1]-v[i-1]) / dt;
    }

    // Upstroke: first index where dV/dt crosses dvdt_threshold upward
    size_t i_up = n; // invalid
    for (size_t i=1; i<n; ++i)
    {
        if (dvdt[i-1] < dvdt_threshold && dvdt[i] >= dvdt_threshold) { i_up = i; break; }
    }
    if (i_up == n) return bm;
    bm.t_up_ms = t[i_up];

    // Diastolic window: 20 ms before upstroke (clamped to start)
    double t_dia_start = std::max(t.front(), bm.t_up_ms - 20.0);
    double rmp = 0.0; unsigned rmp_count=0;
    for (size_t i=0; i<n; ++i)
    {
        if (t[i] >= t_dia_start && t[i] < bm.t_up_ms) { rmp += v[i]; ++rmp_count; }
    }
    bm.rmp_mV = (rmp_count ? rmp / rmp_count : v.front());

    // Peak and APA
    auto it_peak = std::max_element(v.begin()+i_up, v.end());
    size_t i_peak = static_cast<size_t>(std::distance(v.begin(), it_peak));
    bm.vpeak_mV = *it_peak;
    bm.apa_mV   = bm.vpeak_mV - bm.rmp_mV;

    // dV/dt max around upstroke (pm5 ms)
    const double t_lo = bm.t_up_ms - 5.0, t_hi = bm.t_up_ms + 5.0;
    double dvdt_max = dvdt[i_up];
    for (size_t i=1; i+1<n; ++i)
    {
        if (t[i] >= t_lo && t[i] <= t_hi) dvdt_max = std::max(dvdt_max, dvdt[i]);
    }
    bm.dvdt_max = dvdt_max;

    // APD levels
    const double v50 = bm.vpeak_mV - 0.5*bm.apa_mV;
    const double v90 = bm.vpeak_mV - 0.9*bm.apa_mV;

    // Find repolarization crossings after peak
    auto find_cross = [&](double vlevel)->double {
        for (size_t i=i_peak+1; i<n; ++i)
        {
            if (v[i-1] >= vlevel && v[i] < vlevel)
            {
                // linear interpolate
                const double a = v[i-1]-vlevel, b = v[i]-vlevel;
                const double frac = a / (a - b + 1e-12);
                return t[i-1] + frac*(t[i]-t[i-1]);
            }
        }
        return std::numeric_limits<double>::quiet_NaN();
    };

    const double t50 = find_cross(v50);
    const double t90 = find_cross(v90);

    if (!std::isnan(t50)) bm.apd50_ms = t50 - bm.t_up_ms;
    if (!std::isnan(t90)) bm.apd90_ms = t90 - bm.t_up_ms;
    if (!std::isnan(bm.apd90_ms) && !std::isnan(bm.apd50_ms))
        bm.triangulation_ms = bm.apd90_ms - bm.apd50_ms;

    return bm;
}

/*
        Simulate one beat and sample Vm uniformly with dt_sample.
        Returns (t, v) over [t0, t0+period].
 */
inline void SimulateOneBeat(AbstractCvodeCell& cell,
                            double t0_ms, double period_ms, double dt_sample_ms,
                            std::vector<double>& t, std::vector<double>& v)
{
    t.clear(); v.clear();
    const unsigned ns = static_cast<unsigned>(std::ceil(period_ms / dt_sample_ms));
    double tcur = t0_ms;
    double tend = t0_ms + period_ms;
    t.reserve(ns+1); v.reserve(ns+1);

    // Ensure CVODE cell starts from t0 (if needed)
    cell.SetTimestep(dt_sample_ms);
    // Sample at uniform grid
    for (;;)
    {
        t.push_back(tcur);
        v.push_back(cell.GetVoltage());
        if (tcur >= tend) break;
        const double tnext = std::min(tcur + dt_sample_ms, tend);
        cell.SolveAndUpdateState(tcur, tnext);
        tcur = tnext;
    }
}

/*
        Pace to steady state and compute biomarkers on the final beat.
 */
inline void RunSingleCell(const SingleCellConfig& cfg)
{
    OutputFileHandler out(cfg.outdir, /*cleanOutputDirectory=*/true);

    // Build model + pacing
    boost::shared_ptr<RegularStimulus> stim(new RegularStimulus(
        cfg.stim_amp_uAcm2, cfg.stim_dur_ms, cfg.period_ms, cfg.stim_start_ms));

    pig_ventr_ap_endoCvodeOpt cell(boost::shared_ptr<AbstractIvpOdeSolver>(), stim);

    // Beat loop
    std::vector<double> apd90_hist;
    apd90_hist.reserve(cfg.max_beats);

    std::vector<double> t_last, v_last;
    BeatMetrics bm_last;

    double t0 = 0.0;
    for (unsigned beat=1; beat<=cfg.max_beats; ++beat)
    {
        // Simulate one full period and sample Vm
        SimulateOneBeat(cell, t0, cfg.period_ms, cfg.dt_sample_ms, t_last, v_last);
        t0 += cfg.period_ms;

        // Analyze this beat's APD90 (for convergence check)
        BeatMetrics bm = AnalyzeBeat(t_last, v_last, cfg.dvdt_threshold);
        apd90_hist.push_back(bm.apd90_ms);
        bm_last = bm; // keep last

        // Check steady state after minimum number of beats
        if (beat >= std::max(cfg.min_beats, cfg.steady_window))
        {
            const unsigned K = cfg.steady_window;
            bool ok = true;
            // Require all pairwise diffs within last K beats to be <= tol (robust)
            for (unsigned i=0; i<K && ok; ++i)
            {
                for (unsigned j=i+1; j<K && ok; ++j)
                {
                    const double a = apd90_hist[apd90_hist.size()-1-i];
                    const double b = apd90_hist[apd90_hist.size()-1-j];
                    if (std::isnan(a) || std::isnan(b) || std::fabs(a-b) > cfg.apd_tol_ms) ok = false;
                }
            }
            if (ok) break;
        }
    }

    // Final biomarkers on the last beat (already computed as bm_last)
    // Also compute RMP more robustly from late-diastolic window of the last beat
    // (AnalyzeBeat already averages 20 ms before upstroke).

    // Write biomarkers CSV
    {
        std::ostringstream ss;
        ss << "period_ms, stim_start_ms, stim_dur_ms, stim_amp_uAcm2,"
              "RMP_mV, Vpeak_mV, APA_mV, dVdt_max_mV_per_ms, APD50_ms, APD90_ms, Triang_ms, t_up_ms\n";
        ss << std::fixed << std::setprecision(6)
           << cfg.period_ms << "," << cfg.stim_start_ms << "," << cfg.stim_dur_ms << ","
           << cfg.stim_amp_uAcm2 << ","
           << bm_last.rmp_mV << "," << bm_last.vpeak_mV << "," << bm_last.apa_mV << ","
           << bm_last.dvdt_max << "," << bm_last.apd50_ms << "," << bm_last.apd90_ms << ","
           << bm_last.triangulation_ms << "," << bm_last.t_up_ms << "\n";
        WriteCsv(out.GetOutputDirectoryFullPath() + "singlecell_biomarkers.csv", ss.str());
    }

    // Write final beat trace
    if (cfg.write_trace_csv)
    {
        std::ostringstream ts;
        ts << "time_ms,Vm_mV\n";
        for (size_t i=0; i<t_last.size(); ++i)
            ts << std::fixed << std::setprecision(6) << t_last[i] << "," << v_last[i] << "\n";
        WriteCsv(out.GetOutputDirectoryFullPath() + "singlecell_trace.csv", ts.str());
    }

    if (PetscTools::AmMaster())
    {
        std::cout << "Wrote biomarkers and trace to: "
                  << out.GetOutputDirectoryFullPath() << std::endl;
    }
}

// -------- main --------
#ifdef PORCISA_SINGLECELL_MAIN
int main(int argc, char* argv[])
{
    CommandLineArguments* p_args = CommandLineArguments::Instance();
    p_args->Parse(argc, argv);

    SingleCellConfig cfg;

    if (p_args->OptionExists("--out"))     cfg.outdir = p_args->GetOption("--out");
    if (p_args->OptionExists("--bcl"))     cfg.period_ms = atof(p_args->GetOption("--bcl").c_str());
    if (p_args->OptionExists("--stimt"))   cfg.stim_start_ms = atof(p_args->GetOption("--stimt").c_str());
    if (p_args->OptionExists("--stimdur")) cfg.stim_dur_ms = atof(p_args->GetOption("--stimdur").c_str());
    if (p_args->OptionExists("--stimamp")) cfg.stim_amp_uAcm2 = atof(p_args->GetOption("--stimamp").c_str());
    if (p_args->OptionExists("--maxb"))    cfg.max_beats = (unsigned)atoi(p_args->GetOption("--maxb").c_str());
    if (p_args->OptionExists("--minb"))    cfg.min_beats = (unsigned)atoi(p_args->GetOption("--minb").c_str());
    if (p_args->OptionExists("--kwin"))    cfg.steady_window = (unsigned)atoi(p_args->GetOption("--kwin").c_str());
    if (p_args->OptionExists("--apdtol"))  cfg.apd_tol_ms = atof(p_args->GetOption("--apdtol").c_str());
    if (p_args->OptionExists("--dt"))      cfg.dt_sample_ms = atof(p_args->GetOption("--dt").c_str());
    if (p_args->OptionExists("--dvdt"))    cfg.dvdt_threshold = atof(p_args->GetOption("--dvdt").c_str());
    if (p_args->OptionExists("--trace"))   cfg.write_trace_csv = (atoi(p_args->GetOption("--trace").c_str()) != 0);

    try {
        RunSingleCell(cfg);
    } catch (const std::exception& e) {
        if (PetscTools::AmMaster())
            std::cerr << "Fatal: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
#endif // PORCISA_SINGLECELL_MAIN

#endif // TEST_SINGLE_CELL_HPP_
