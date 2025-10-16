#include <fstream>

#include "SimpleStimulus.hpp"

#include "carro_2011_epi_pigCvodeOpt.hpp"
#include "carro_2011_endo_pigCvodeOpt.hpp"
#include "carro_2011_m_pigCvodeOpt.hpp"

class PorcineCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    PorcineHeterogen* mpPorcineHeterogen;
    std::vector<boost::shared_ptr<SimpleStimulus> > mSimpleStimuli;
    std::vector<double> mStimulusinfo;
	bool mIschaemia;
	
	double	_k0;
	double	_ikatp;
	double	_ina;
	double	_ica;
	double  _gncx;
    double  _pnak;
    double  _inal;
	double	_radius;
	double	_xx;
	double	_yy;
	double	_zz;	
	
	unsigned _transmural;														

	double mLayer1, mLayer2;
	unsigned mOptionStimulus;

public:
    PorcineCellFactory(
    	const std::vector<double>& stimulusinfo,
		std::string meshpath,
		PorcineHeterogen* pPorcineHeterogen,
		const std::vector<unsigned>& rNodePerm,
		bool ischaemia,
		double	k0,
		double	ikatp,
		double	ina,
		double	ica,
		double  gncx,
        double  pnak,
        double 	inal,
		double	radius,
		double	xx,
		double	yy,
		double	zz,	
		unsigned transmural
		)
							
    : AbstractCardiacCellFactory<3>(),
    
    mpPorcineHeterogen(pPorcineHeterogen),
    mStimulusinfo(stimulusinfo),
    _k0(k0),
	_ikatp(ikatp),
	_ina(ina),
	_ica(ica),
	_gncx(gncx),
    _pnak(pnak),
    _inal(inal),
	_radius(radius),
	_xx(xx),
	_yy(yy),
	_zz(zz),
	_transmural(transmural)

	{
		mIschaemia = ischaemia;	
	
		mLayer1 = double(cellcomposition[0]) / 100.0f;
		mLayer2 = double(cellcomposition[1]) / 100.0f;
		
		mOptionStimulus = mStimulusinfo.size();
	}

    AbstractCardiacCellInterface* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        AbstractCvodeCell* p_cell;

        unsigned node_index = pNode->GetIndex();

        HeartGeometryInformation<3>* p_heart_geom_info = GetHeartGeometryInformation();
        const double lv_dist = p_heart_geom_info->rGetDistanceMapLeftVentricle()[node_index];
        const double rv_dist = p_heart_geom_info->rGetDistanceMapRightVentricle()[node_index];

        double distance_epi = p_heart_geom_info->rGetDistanceMapEpicardium()[node_index];
        double distance_endo = std::min(lv_dist, rv_dist);
        double relative_position = distance_endo / (distance_endo + distance_epi);
		
		if (relative_position > (1.0 - mLayer1))
			{
				p_cell = new carro_2011_epi_pigCvodeOpt(mpSolver, mpZeroStimulus);					
			}
		else if(relative_position > (1.0 - mLayer1 - mLayer2))
			{
				p_cell = new carro_2011_m_pigCvodeOpt(mpSolver, mpZeroStimulus);
			}
		else
			{
				p_cell = new carro_2011_endo_pigCvodeOpt(mpSolver, mpZeroStimulus);
			}
		
		double isch_r = _radius;
        double isch_bz_r = isch_r + 0.5;
        double x_isch = _xx;double y_isch = _yy; double z_isch=_zz;
        ChastePoint<3> isch_centre (x_isch,y_isch,z_isch);
        ChastePoint<3> isch_radius (isch_r,isch_r,isch_r);
        ChastePoint<3> isch_bz_radius (isch_bz_r,isch_bz_r,isch_bz_r);
        ChasteEllipsoid<3> isch_region (isch_centre,isch_radius);
        ChasteEllipsoid<3> isch_bz_region (isch_centre,isch_bz_radius);

        bool cell_is_in_ischReg = isch_region.DoesContain(this->GetMesh()->GetNode(node_index)->rGetLocation());
        bool cell_is_in_isch_BZReg = isch_bz_region.DoesContain(this->GetMesh()->GetNode(node_index)->rGetLocation());
	
        double x = this->GetMesh()->GetNode(node_index)->rGetLocation()[0];
        double y = this->GetMesh()->GetNode(node_index)->rGetLocation()[1];
        double z = this->GetMesh()->GetNode(node_index)->rGetLocation()[2];
		
		// ischaemic parameters
				
		double k_n = 5;
		double k_i = _k0;
		double fkatp_n = 0;
		double fkatp_i = _ikatp;
		double ina_n = 1.0;
		double ina_i = _ina;
		double ica_n = 1.0;
		double ica_i = _ica;

		double gncx_n = 1;
        double gncx_i = _gncx;
        double pnak_n = 1;
        double pnak_i = _pnak;

        double inal_n = 1;
        double inal_i = _inal;
		
		bool limit;
		
		if (_transmural == 1)
			limit = relative_position > 0.1;
			
		if (_transmural == 2)
			limit = (relative_position > 0.1) && (relative_position < 0.3);
			
		if (_transmural == 0)
			limit = (relative_position > 0.1) && (relative_position < 0.5);
		
		if (cell_is_in_isch_BZReg && mIschaemia) 
		{

			if(limit)
			{
		
				if(cell_is_in_ischReg) // ICZ
				{
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

				else if (_transmural==1) // BZ
				{
					double distance_from_isch_centre = sqrt(pow((x-x_isch),2)+pow((y-y_isch),2)+pow((z-z_isch),2));
					double distance_from_cz = distance_from_isch_centre - isch_r;

					double bz_k = isch_bz_r-isch_r; // K+ BZ width

					double bz_ina = 0.5*bz_k; // INa and ICaL BZ width
					double distance_from_cz_ina = distance_from_cz - bz_ina;
					
					double bz_fkatp = bz_ina;
					double distance_from_cz_fkatp = distance_from_cz - bz_fkatp;

					double bz_gncx = bz_ina; // NCX BZ width
                    double distance_from_cz_gncx = distance_from_cz - bz_gncx;
                    double bz_pnak = bz_ina; //bz_fkatp; // pnak BZ width
                    double distance_from_cz_pnak = distance_from_cz - bz_pnak;

                    double bz_inal = bz_ina;
					double distance_from_cz_inal = distance_from_cz - bz_inal;

					double a_k = (k_n-k_i)/bz_k;
					double a_ina = (ina_n-ina_i)/bz_ina;
					double a_ica = (ica_n-ica_i)/bz_ina;
					double a_fkatp = (fkatp_n-fkatp_i)/bz_fkatp;

					double a_gncx = (gncx_n-gncx_i)/bz_gncx;
                    double a_pnak = (pnak_n-pnak_i)/bz_pnak;

                    double a_inal = (inal_n-inal_i)/bz_inal;

					double y_k = a_k*distance_from_cz+k_i;
					double y_ina = a_ina*distance_from_cz_ina+ina_i;
					double y_ica = a_ica*distance_from_cz_ina+ica_i;
					double y_fkatp = a_fkatp*distance_from_cz_fkatp+fkatp_i;

					double y_gncx = a_gncx*distance_from_cz_gncx+gncx_i;
                    double y_pnak = a_pnak*distance_from_cz_pnak+pnak_i;

                    double y_inal = a_inal*distance_from_cz_inal+inal_i;

					p_cell->SetParameter("extracellular_potassium_concentration",y_k);
					double gna = p_cell->GetParameter("membrane_fast_sodium_current_conductance");
					double gca = p_cell->GetParameter("membrane_L_type_calcium_current_conductance");
					p_cell->SetParameter("membrane_fast_sodium_current_conductance",gna*y_ina);
					p_cell->SetParameter("membrane_L_type_calcium_current_conductance",gca*y_ica);
					p_cell->SetParameter("membrane_atp_dependent_potassium_current_conductance",y_fkatp);

					double get_gNCX = p_cell->GetParameter("membrane_sodium_calcium_exchanger_current_conductance");
                    p_cell->SetParameter("membrane_sodium_calcium_exchanger_current_conductance", get_gNCX * y_gncx);
                    double get_pNaK = p_cell->GetParameter("membrane_sodium_potassium_pump_current_permeability");
                    p_cell->SetParameter("membrane_sodium_potassium_pump_current_permeability", get_pNaK * y_pnak);

                    double get_gNaL = p_cell->GetParameter("membrane_late_sodium_current_conductance");
                    p_cell->SetParameter("membrane_late_sodium_current_conductance", get_gNaL * y_inal);

				}
			}
		}
				
		// applying stimuli

		double extra_small_z = 0.25;

        double stim_amplitude = -45000;
		double stim_duration = 2;

		double  ax= 21.8 - 0.264 ,   ay= 6.9 + 0.793,   az= 54.5 - 0.83 - 0.3; // basal septum (5 to 17 ms)
		double  bx= 21.8 - 0.21 ,   by= 6.9 - 0.59,   bz= 54.5  - 2.46 - 0.3; // central septum (7-10 ms)
		double  cx= 21.8 +  1.36 ,   cy= 6.9 - 2.55  ,   cz= 54.5  -  4.9 - 0.2; // anterior apex(11-15-18 ms)
		double  dx= 21.8 + 1.4,   dy= 6.9 +2.07,   dz= 54.5 -3.35 - 0.3; //  posterior central (approx 23 ms)
		double  ex= 21.8 + 1.84,   ey= 6.9 + 1.09,   ez= 54.5 - 5.25 - 0.2 ; //  posterior apex (15-20ms)
		double  fx= 21.8 + 2.25,   fy= 6.9 - 2.16,   fz= 54.5 - 2.6 - 0.3; //  anterior central (approx 15 ms) 
		double  gx= 21.8 + 3.6,   gy= 6.9 + 0.66,   gz= 54.5 - 1.58 - 0.3; //  lv free wall central (approx 22 ms)
		double  hx= 21.8 - 1.7,   hy= 6.9 - 2.38,   hz= 54.5 - 2.626 - 0.3; //  rv posterior apex (approx 11-15 ms) 	

		ChastePoint<3>  a_root_node_c   (   ax  ,   ay  ,   az  )   ;
		ChastePoint<3>  b_root_node_c   (   bx  ,   by  ,   bz  )   ;
		ChastePoint<3>  c_root_node_c   (   cx  ,   cy  ,   cz  )   ;
		ChastePoint<3>  d_root_node_c   (   dx  ,   dy  ,   dz  )   ;
		ChastePoint<3>  e_root_node_c   (   ex  ,   ey  ,   ez  )   ;
		ChastePoint<3>  f_root_node_c   (   fx  ,   fy  ,   fz  )   ;
		ChastePoint<3>  g_root_node_c   (   gx  ,   gy  ,   gz  )   ;
		ChastePoint<3>  h_root_node_c   (   hx  ,   hy  ,   hz  )   ;
		                                   
		ChastePoint<3>  a_root_node_r   (   extra_small_z  +0.1  ,   extra_small_z  +0.1  ,   extra_small_z  +0.1  )   ;
		ChastePoint<3>  b_root_node_r   (   extra_small_z  +0.1 ,   extra_small_z  +0.1 ,   extra_small_z  +0.1 )   ;
		ChastePoint<3>  c_root_node_r   (   extra_small_z  +0.1 ,   extra_small_z  +0.1 ,   extra_small_z  +0.1 )   ;
		ChastePoint<3>  d_root_node_r   (   extra_small_z   +0.1,   extra_small_z   +0.1,   extra_small_z   +0.1)   ;
		ChastePoint<3>  e_root_node_r   (   extra_small_z  +0.1 ,   extra_small_z   +0.1,   extra_small_z  +0.1 )   ;
		ChastePoint<3>  f_root_node_r   (   extra_small_z   +0.1,   extra_small_z   +0.1,   extra_small_z  +0.1 )   ;
		ChastePoint<3>  g_root_node_r   (   extra_small_z  +0.1 ,   extra_small_z   +0.1,   extra_small_z  +0.1 )   ;
		ChastePoint<3>  h_root_node_r   (   extra_small_z  +0.1 ,   extra_small_z   +0.1,   extra_small_z  +0.1 )   ;
		                            
		ChasteEllipsoid<3>  a_rut   (   a_root_node_c   ,   a_root_node_r   )   ;
		ChasteEllipsoid<3>  b_rut   (   b_root_node_c   ,   b_root_node_r   )   ;
		ChasteEllipsoid<3>  c_rut   (   c_root_node_c   ,   c_root_node_r   )   ;
		ChasteEllipsoid<3>  d_rut   (   d_root_node_c   ,   d_root_node_r   )   ;
		ChasteEllipsoid<3>  e_rut   (   e_root_node_c   ,   e_root_node_r   )   ;
		ChasteEllipsoid<3>  f_rut   (   f_root_node_c   ,   f_root_node_r   )   ;
		ChasteEllipsoid<3>  g_rut   (   g_root_node_c   ,   g_root_node_r   )   ;
		ChasteEllipsoid<3>  h_rut   (   h_root_node_c   ,   h_root_node_r   )   ;
		               
		bool    a_rootnode  =   a_rut.DoesContain(this->GetMesh()->GetNode(node_index)->rGetLocation());
		bool    b_rootnode  =   b_rut.DoesContain(this->GetMesh()->GetNode(node_index)->rGetLocation());
		bool    c_rootnode  =   c_rut.DoesContain(this->GetMesh()->GetNode(node_index)->rGetLocation());
		bool    d_rootnode  =   d_rut.DoesContain(this->GetMesh()->GetNode(node_index)->rGetLocation());
		bool    e_rootnode  =   e_rut.DoesContain(this->GetMesh()->GetNode(node_index)->rGetLocation());
		bool    f_rootnode  =   f_rut.DoesContain(this->GetMesh()->GetNode(node_index)->rGetLocation());
		bool    g_rootnode  =   g_rut.DoesContain(this->GetMesh()->GetNode(node_index)->rGetLocation());
		bool    h_rootnode  =   h_rut.DoesContain(this->GetMesh()->GetNode(node_index)->rGetLocation());

		boost::shared_ptr<SimpleStimulus> p_stim_a(new SimpleStimulus(stim_amplitude, stim_duration, mStimulusinfo[0])); //
		boost::shared_ptr<SimpleStimulus> p_stim_b(new SimpleStimulus(stim_amplitude, stim_duration, mStimulusinfo[0] + 3.5)); // + 3.5
		boost::shared_ptr<SimpleStimulus> p_stim_c(new SimpleStimulus(stim_amplitude, stim_duration, mStimulusinfo[0] + 6 )); // + 6 
		boost::shared_ptr<SimpleStimulus> p_stim_d(new SimpleStimulus(stim_amplitude, stim_duration, mStimulusinfo[0] + 16)); // + 16
		boost::shared_ptr<SimpleStimulus> p_stim_e(new SimpleStimulus(stim_amplitude, stim_duration, mStimulusinfo[0] +8 )); //  +8
		boost::shared_ptr<SimpleStimulus> p_stim_f(new SimpleStimulus(stim_amplitude, stim_duration, mStimulusinfo[0] +18)); // + 18
		boost::shared_ptr<SimpleStimulus> p_stim_g(new SimpleStimulus(stim_amplitude, stim_duration, mStimulusinfo[0] + 20)); // + 20
		boost::shared_ptr<SimpleStimulus> p_stim_h(new SimpleStimulus(stim_amplitude, stim_duration, mStimulusinfo[0] + 22)); //+ 22

                                if ( (a_rootnode) )
                                {
                                    p_cell->SetStimulusFunction(p_stim_a);
                                }
                                if ( (b_rootnode) )
                                {
                                    p_cell->SetStimulusFunction(p_stim_b);
                                }
                                if ( (c_rootnode) )
                                {
                                    p_cell->SetStimulusFunction(p_stim_c);
                                }
                                if ( (d_rootnode) )
                                {
                                    p_cell->SetStimulusFunction(p_stim_d);
                                }
                                if ( (e_rootnode) )
                                {
                                    p_cell->SetStimulusFunction(p_stim_e);
                                }
                                if ( (f_rootnode) )
                                {
                                    p_cell->SetStimulusFunction(p_stim_f);
                                }
                                if ( (g_rootnode) )
                                {
                                    p_cell->SetStimulusFunction(p_stim_g);
                                }
                                if ( (h_rootnode) )
                                {
                                    p_cell->SetStimulusFunction(p_stim_h);
                                }

 
		const double gks = p_cell->GetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance");
        const double gks_scaling = mpPorcineHeterogen->GetScalingFactor(*pNode);
        p_cell->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", gks_scaling*gks);
		p_cell->SetTolerances(1e-4,1e-6);

		return p_cell;
    }

};
