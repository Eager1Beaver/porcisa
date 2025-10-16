#include <cxxtest/TestSuite.h>
#include <ctime>
#include "PetscSetupAndFinalize.hpp"

#include "HeartEventHandler.hpp"

#include "SingleTraceOutputModifier.hpp"
#include "ActivationOutputModifier.hpp" 

#include "BidomainWithBathProblem.hpp"
#include "HeartGeometryInformation.hpp"

#include "PorcineHeterogen.hpp"
#include "PorcineCellFactory.hpp"

#include "PorcineConductivityModifier.hpp"

class TestPorcineECG : public CxxTest::TestSuite
{

public:
    void TestMain()
    {
        HeartEventHandler::BeginEvent(HeartEventHandler::EVERYTHING);
        
        ///////////////////////////////////////////////////////////////////////////////////////     
        // Parsing arguments //////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////// 

        	double odet = 0.1; 	
        	double pdet = 0.1;
        	double simduration = 350; 
			
			std::string meshname 		= 	"pig91";

			std::string preconditioner;
			std::string solver;

			solver = "symmlq"; 

			bool regularstim = false, simplestim = false; 
			bool ischaemia = false, phase1a = false, phase1b = false, phase1a1b = false, phase0a = false;
			bool ectopic = false;
			double period;
			double simplestimtime;
			//double ectopicS2;

			bool post_process = true; 	
			bool getAPD = true;

			bool ischaemia_phase1a = false, ischaemia_phase1b = false, ischaemia_phase1a1b = false, ischaemia_phase0a = false;

			std::string phase = "ctrl";

			if (ischaemia_phase1a == true)
			{
				phase = "ph1a";
				ischaemia = true;
				phase1a = true;
			}

			if (ischaemia_phase1b == true)
			{
				phase = "ph1b";
				ischaemia = true;
				phase1b = true;
			}

			if (ischaemia_phase1a1b == true)
			{
				phase = "ph1a1b";
				ischaemia = true;
				phase1a1b = true;
			}

			if (ischaemia_phase0a == true)
			{
				phase = "ph0a";
				ischaemia = true;
				phase0a = true;
			}

			std::string output_directory = 	meshname + "_" + phase + "_N_ctrl";		
        	std::string filename_prefix = 	"pecg_" + phase;					
		
			double k0		= 0.0;
			double ikatp	= 0.0;
			double ina		= 0.0;
			double ica		= 0.0;
			double gncx     = 0.0; 	
        	double pnak     = 0.0; 		
        	double inal		= 0.0;		
        	double factor_h	= 1.0;	
			double radius	= 0.0;
			double xx		= 0.0;
			double yy		= 0.0;
			double zz		= 0.0;
			
			//double s2radius  = 0.0;
			//double s2x		 = 0.0;
			//double s2y		 = 0.0;
			//double s2z		 = 0.0;
			unsigned transmural = 1;
			double printsteps = 1.0; 	
						

			regularstim=true; 	
			period = 1000; 		

			if (ischaemia == true)
			{
				if (phase1a == true)
				{
					k0 = 	6.5;
					ikatp = 0.06; 
					ina = 	0.85;
					ica = 	0.85; 
					inal = 	1.6; 
					gncx = 1;
					pnak = 1;

					factor_h = 0.22;
				}

				if (phase1b == true) 
				{
					k0 = 	8.3; 
					ikatp = 0.1; 
					ina = 	0.5; 
					ica = 	0.5; 
					inal = 	3.0;
					gncx = 	0.33;
					pnak = 	0.33; 

					factor_h = 0.1;
				}

				if (phase1a1b == true)
				{
					k0 = 	8.1; 
					ikatp = 0.08; 
					ina = 	0.65;
					ica = 	0.65; 
					inal = 2.8; 
					gncx = 0.5; 
					pnak = 0.5; 

					factor_h = 0.18; 

				if (phase0a == true)
				{
					k0 = 	5.5; 
					ikatp = 0.01; 
					ina = 	0.92; 
					ica = 	0.92; 
					inal = 	1.2; 
					gncx = 1;
					pnak = 1;

					factor_h = 0.5;
				}

			}

				radius			= 2.3; 
				xx 				= 23;
				yy 				= 4.9; 
				zz 				= 49.15;	
				transmural 		= 1;

				/*if (ectopicS2)
				{
					ectopic=true;
					
					s2radius	= 0;
					s2x			= 0;
					s2y			= 0;
					s2z			= 0;
				}*/
				
				preconditioner = "blockdiagonal"; 

			unsigned epi_layer 			= 25; 		
			unsigned endo_layer 		= 30; 			

		///////////////////////////////////////////////////////////////////////////////////////     
        // Mesh ///////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////// 

			std::string filepath = "./meshes/" + meshname + "/" + meshname + "_bin";
			HeartConfig::Instance()->SetMeshFileName(filepath, cp::media_type::Orthotropic); 

			// terminal outputting
			//std::cout << "loading mesh from: " << filepath << std::endl;
			std::cout << "odet: " << odet << ", pdet: " << pdet << ", simduration: " << simduration << std::endl;
			std::cout << "solver: " << solver << ", preconditioner: " << preconditioner << std::endl;
			std::cout << "ischaemia: " << ischaemia << ", phase1a: " << phase1a << ", phase1b: " << phase1b << ", phase1a1b: " << phase1a1b << ", phase0a: " << phase0a << std::endl;
			std::cout << "factor_h: " << factor_h << std::endl;
			std::cout << "output_directory: " << output_directory << std::endl;
						
			// Reading mesh
			HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
			DistributedTetrahedralMesh<3,3> mesh;
			TrianglesMeshReader<3,3> mesh_reader(filepath);
			mesh.ConstructFromMeshReader(mesh_reader);
			HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);
			HeartEventHandler::BeginEvent(HeartEventHandler::INITIALISE);

			// Heart geometry information
			const std::string epi_face_file = filepath+".epi";
			const std::string rv_face_file = filepath+".rv";
			const std::string lv_face_file = filepath+".lv";
			HeartGeometryInformation<3> heart_geom_info(mesh, epi_face_file, lv_face_file, rv_face_file, true);

        ///////////////////////////////////////////////////////////////////////////////////////     
        // Output /////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////
			
			HeartConfig::Instance()->SetVisualizeWithMeshalyzer(false);
			HeartConfig::Instance()->SetVisualizeWithCmgui(false);
			HeartConfig::Instance()->SetVisualizeWithVtk(false);
			HeartConfig::Instance()->SetVisualizeWithParallelVtk(true);       
			
        ///////////////////////////////////////////////////////////////////////////////////////     
        // Numerical setup ////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////// 
        		
			HeartConfig::Instance()->SetUseStateVariableInterpolation(false);

			HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(odet, pdet, printsteps);
			
			HeartConfig::Instance()->SetKSPPreconditioner( preconditioner.c_str() );
			HeartConfig::Instance()->SetKSPSolver( solver.c_str() );

        ///////////////////////////////////////////////////////////////////////////////////////     
        // Simulation setup ///////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////// 
        
        
			// Simulation duration
			HeartConfig::Instance()->SetSimulationDuration(simduration); // ms
			
			// Conductivities

			double intra_longi = 0 , intra_trans = 0, intra_normal = 0;
        	double extra_longi = 0 , extra_trans = 0, extra_normal = 0;

        	intra_longi = 3.35; intra_trans = 1.755; intra_normal = 0.8775;
        	extra_longi = 12.4; extra_trans = 5.5; extra_normal = 5.5;

        	HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(intra_longi, intra_trans, intra_normal));
        	HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(extra_longi, extra_trans, extra_normal));
        		
        		if (ischaemia == true)
        		{
        			double isch_radius_h = radius;

        			std::vector<ChasteEllipsoid<3> > heterogeneity_isch_area0;
        			std::vector<ChasteEllipsoid<3> > heterogeneity_isch_area1;
        			std::vector<ChasteEllipsoid<3> > heterogeneity_isch_area2;
        			std::vector<ChasteEllipsoid<3> > heterogeneity_isch_area3;
        			std::vector<ChasteEllipsoid<3> > heterogeneity_isch_area4;
        			std::vector<ChasteEllipsoid<3> > heterogeneity_isch_area5;
        			std::vector<ChasteEllipsoid<3> > heterogeneity_isch_area6;
        			std::vector<ChasteEllipsoid<3> > heterogeneity_isch_area7;
        			std::vector<ChasteEllipsoid<3> > heterogeneity_isch_area8;
        			std::vector<ChasteEllipsoid<3> > heterogeneity_isch_area9;
        			std::vector<ChasteEllipsoid<3> > heterogeneity_isch_area10;

        			std::vector< c_vector<double,3> > intra_isch_conductivities0;
        			std::vector< c_vector<double,3> > intra_isch_conductivities1;
        			std::vector< c_vector<double,3> > intra_isch_conductivities2;
        			std::vector< c_vector<double,3> > intra_isch_conductivities3;
        			std::vector< c_vector<double,3> > intra_isch_conductivities4;
        			std::vector< c_vector<double,3> > intra_isch_conductivities5;
        			std::vector< c_vector<double,3> > intra_isch_conductivities6;
        			std::vector< c_vector<double,3> > intra_isch_conductivities7;
        			std::vector< c_vector<double,3> > intra_isch_conductivities8;
        			std::vector< c_vector<double,3> > intra_isch_conductivities9;
        			std::vector< c_vector<double,3> > intra_isch_conductivities10;

        			ChastePoint<3> isch_centre_h(xx, yy, zz); // common

        			ChastePoint<3> isch_radii_h0(isch_radius_h, isch_radius_h, isch_radius_h);
        			ChastePoint<3> isch_radii_h1(isch_radius_h+0.05, isch_radius_h+0.05, isch_radius_h+0.05);
        			ChastePoint<3> isch_radii_h2(isch_radius_h+0.1, isch_radius_h+0.1, isch_radius_h+0.1);
        			ChastePoint<3> isch_radii_h3(isch_radius_h+0.15, isch_radius_h+0.15, isch_radius_h+0.15);
        			ChastePoint<3> isch_radii_h4(isch_radius_h+0.2, isch_radius_h+0.2, isch_radius_h+0.2);
        			ChastePoint<3> isch_radii_h5(isch_radius_h+0.25, isch_radius_h+0.25, isch_radius_h+0.25);
        			ChastePoint<3> isch_radii_h6(isch_radius_h+0.3, isch_radius_h+0.3, isch_radius_h+0.3);
        			ChastePoint<3> isch_radii_h7(isch_radius_h+0.35, isch_radius_h+0.35, isch_radius_h+0.35);
        			ChastePoint<3> isch_radii_h8(isch_radius_h+0.40, isch_radius_h+0.40, isch_radius_h+0.40);
        			ChastePoint<3> isch_radii_h9(isch_radius_h+0.45, isch_radius_h+0.45, isch_radius_h+0.45);
        			ChastePoint<3> isch_radii_h10(isch_radius_h+0.5, isch_radius_h+0.5, isch_radius_h+0.5);

        			ChasteEllipsoid<3> isch_region_h0(isch_centre_h, isch_radii_h0);
        			ChasteEllipsoid<3> isch_region_h1(isch_centre_h, isch_radii_h1);
        			ChasteEllipsoid<3> isch_region_h2(isch_centre_h, isch_radii_h2);
        			ChasteEllipsoid<3> isch_region_h3(isch_centre_h, isch_radii_h3);
        			ChasteEllipsoid<3> isch_region_h4(isch_centre_h, isch_radii_h4);
        			ChasteEllipsoid<3> isch_region_h5(isch_centre_h, isch_radii_h5);
        			ChasteEllipsoid<3> isch_region_h6(isch_centre_h, isch_radii_h6);
        			ChasteEllipsoid<3> isch_region_h7(isch_centre_h, isch_radii_h7);
        			ChasteEllipsoid<3> isch_region_h8(isch_centre_h, isch_radii_h8);
        			ChasteEllipsoid<3> isch_region_h9(isch_centre_h, isch_radii_h9);
        			ChasteEllipsoid<3> isch_region_h10(isch_centre_h, isch_radii_h10);

        			heterogeneity_isch_area0.push_back(isch_region_h0);
        			heterogeneity_isch_area1.push_back(isch_region_h1);
        			heterogeneity_isch_area2.push_back(isch_region_h2);
        			heterogeneity_isch_area3.push_back(isch_region_h3);
        			heterogeneity_isch_area4.push_back(isch_region_h4);
        			heterogeneity_isch_area5.push_back(isch_region_h5);
        			heterogeneity_isch_area6.push_back(isch_region_h6);
        			heterogeneity_isch_area7.push_back(isch_region_h7);
        			heterogeneity_isch_area8.push_back(isch_region_h8);
        			heterogeneity_isch_area9.push_back(isch_region_h9);
        			heterogeneity_isch_area10.push_back(isch_region_h10);

        			std::vector< c_vector<double,3> > extra_conductivities;
        			extra_conductivities.push_back( Create_c_vector(extra_longi, extra_trans, extra_normal) ); // stays ctrl

        			intra_isch_conductivities0.push_back( Create_c_vector( intra_longi*factor_h, intra_trans*factor_h, intra_normal*factor_h) ); 
			        intra_isch_conductivities1.push_back( Create_c_vector( intra_longi*(factor_h+(1-factor_h)/11), intra_trans*(factor_h+(1-factor_h)/11), intra_normal*(factor_h+(1-factor_h)/11) ));
			        intra_isch_conductivities2.push_back( Create_c_vector( intra_longi*(factor_h+2*(1-factor_h)/11), intra_trans*(factor_h+2*(1-factor_h)/11), intra_normal*(factor_h+2*(1-factor_h)/11) ));
			        intra_isch_conductivities3.push_back( Create_c_vector( intra_longi*(factor_h+3*(1-factor_h)/11), intra_trans*(factor_h+3*(1-factor_h)/11), intra_normal*(factor_h+3*(1-factor_h)/11) ));
			        intra_isch_conductivities4.push_back( Create_c_vector( intra_longi*(factor_h+4*(1-factor_h)/11), intra_trans*(factor_h+4*(1-factor_h)/11), intra_normal*(factor_h+4*(1-factor_h)/11) ));
			        intra_isch_conductivities5.push_back( Create_c_vector( intra_longi*(factor_h+5*(1-factor_h)/11), intra_trans*(factor_h+5*(1-factor_h)/11), intra_normal*(factor_h+5*(1-factor_h)/11) ));
			        intra_isch_conductivities6.push_back( Create_c_vector( intra_longi*(factor_h+6*(1-factor_h)/11), intra_trans*(factor_h+6*(1-factor_h)/11), intra_normal*(factor_h+6*(1-factor_h)/11) ));
			        intra_isch_conductivities7.push_back( Create_c_vector( intra_longi*(factor_h+7*(1-factor_h)/11), intra_trans*(factor_h+7*(1-factor_h)/11), intra_normal*(factor_h+7*(1-factor_h)/11) ));
			        intra_isch_conductivities8.push_back( Create_c_vector( intra_longi*(factor_h+8*(1-factor_h)/11), intra_trans*(factor_h+8*(1-factor_h)/11), intra_normal*(factor_h+8*(1-factor_h)/11) ));
			        intra_isch_conductivities9.push_back( Create_c_vector( intra_longi*(factor_h+9*(1-factor_h)/11), intra_trans*(factor_h+9*(1-factor_h)/11), intra_normal*(factor_h+9*(1-factor_h)/11) ));
			        intra_isch_conductivities10.push_back( Create_c_vector( intra_longi*(factor_h+10*(1-factor_h)/11), intra_trans*(factor_h+10*(1-factor_h)/11), intra_normal*(factor_h+10*(1-factor_h)/11) ));
					
					
			        HeartConfig::Instance()->SetConductivityHeterogeneitiesEllipsoid(
		            		heterogeneity_isch_area10, 
		            		intra_isch_conductivities10, 
		            		extra_conductivities);
			        HeartConfig::Instance()->SetConductivityHeterogeneitiesEllipsoid(
		            		heterogeneity_isch_area9, 
		            		intra_isch_conductivities9, 
		            		extra_conductivities);
			        HeartConfig::Instance()->SetConductivityHeterogeneitiesEllipsoid(
		            		heterogeneity_isch_area8, 
		            		intra_isch_conductivities8, 
		            		extra_conductivities);
			        HeartConfig::Instance()->SetConductivityHeterogeneitiesEllipsoid(
		            		heterogeneity_isch_area7, 
		            		intra_isch_conductivities7, 
		            		extra_conductivities);
			        HeartConfig::Instance()->SetConductivityHeterogeneitiesEllipsoid(
		            		heterogeneity_isch_area6, 
		            		intra_isch_conductivities6, 
		            		extra_conductivities);
			        HeartConfig::Instance()->SetConductivityHeterogeneitiesEllipsoid(
		            		heterogeneity_isch_area5, 
		            		intra_isch_conductivities5, 
		            		extra_conductivities);
			        HeartConfig::Instance()->SetConductivityHeterogeneitiesEllipsoid(
		            		heterogeneity_isch_area4, 
		            		intra_isch_conductivities4, 
		            		extra_conductivities);
			        HeartConfig::Instance()->SetConductivityHeterogeneitiesEllipsoid(
		            		heterogeneity_isch_area3, 
		            		intra_isch_conductivities3, 
		            		extra_conductivities);
			        HeartConfig::Instance()->SetConductivityHeterogeneitiesEllipsoid(
		            		heterogeneity_isch_area2, 
		            		intra_isch_conductivities2, 
		            		extra_conductivities);
			        HeartConfig::Instance()->SetConductivityHeterogeneitiesEllipsoid(
		            		heterogeneity_isch_area1, 
		            		intra_isch_conductivities1, 
		            		extra_conductivities);
			        HeartConfig::Instance()->SetConductivityHeterogeneitiesEllipsoid(
		            		heterogeneity_isch_area0, 
		            		intra_isch_conductivities0, 
		            		extra_conductivities);
			    }            
				

			// tissue
			std::set<unsigned> tissue_ids;
			tissue_ids.insert(0); // heart 	

			// bath
			std::set<unsigned> bath_ids;
			bath_ids.insert(2); // torso
			bath_ids.insert(3); // ribs
			bath_ids.insert(4); // lungs
				
			HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids,bath_ids);

			std::map<unsigned, double> multiple_bath_conductivities;
			multiple_bath_conductivities[2] = 2.29;  // torso
			multiple_bath_conductivities[3] = 0.2;   // ribs
			multiple_bath_conductivities[4] = 0.389; // lungs
				
			HeartConfig::Instance()->SetBathMultipleConductivities(multiple_bath_conductivities);
        
		///////////////////////////////////////////////////////////////////////////////////////     
        // Post Processing ////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////

	        /*std::vector<std::string> output_variables;	
	        //output_variables.push_back("cytosolic_calcium_concentration");
	        output_variables.push_back("membrane_calcium_activated_chloride_current");
	        output_variables.push_back("membrane_slow_delayed_rectifier_potassium_current");
	        output_variables.push_back("membrane_rapid_delayed_rectifier_potassium_current");
	        output_variables.push_back("membrane_inward_rectifier_potassium_current");

	        output_variables.push_back("membrane_L_type_calcium_current");
	        output_variables.push_back("membrane_late_sodium_current");
	        output_variables.push_back("membrane_sodium_calcium_exchanger_current");

	        output_variables.push_back("membrane_atp_dependent_potassium_current");
	        HeartConfig::Instance()->SetOutputVariables(output_variables);*/
	            
	        // Body surface traces
			std::vector<unsigned> electrodes;
			std::string surface_filepath = "/home/user/meshes/" + meshname + "/surface_nodes.txt";

			std::ifstream surface_nodes( surface_filepath.c_str() );
			std::string node; 
			
			while ( surface_nodes >> node )
			{
				electrodes.push_back(std::atoi(node.c_str()));
			}
			if( !surface_nodes.eof() )
			{
				EXCEPTION("Couldn't read surface node file.");
			}
			HeartConfig::Instance()->SetRequestedNodalTimeTraces(electrodes);

			if (post_process == true)
			{	

	        // APD map
				if (getAPD == true)
				{
					std::vector<std::pair<double,double> > apd_map;
					apd_map.push_back(std::pair<double, double>(30.0, 0.0));
					apd_map.push_back(std::pair<double, double>(50.0, 0.0));
					apd_map.push_back(std::pair<double, double>(90.0, 0.0)); 
					HeartConfig::Instance()->SetApdMaps(apd_map);
				}
				
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
				double  ax= 21.8 - 0.264 ,   ay= 6.9 + 0.793,   az= 54.5 - 0.83 - 0.3;
				ChastePoint<3> cv_point(ax,ay,az);
				int cv_point_idx;
				cv_point_idx = mesh.GetNearestNodeIndex(cv_point);
				cv_map.push_back(cv_point_idx);						
				HeartConfig::Instance()->SetConductionVelocityMaps(cv_map);

			}	
	                        
	        Warnings::QuietDestroy(); 
        
        
		///////////////////////////////////////////////////////////////////////////////////////     
        // Problem generation /////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////        

			PorcineHeterogen(&heart_geom_info, true, true, true); 	// heart_geom_info, apicobasal, transmural, interventricular
		
			// Stimulus
				std::vector<double> stimulusinfo;
				
				if (simplestim)
				{
					stimulusinfo.push_back(simplestimtime);
				}
					
				if (regularstim)
				{
					stimulusinfo.push_back(0.0);
					stimulusinfo.push_back(period);
				}
					
				// Cell composition
					std::vector<unsigned> myocard_layer;
				
					myocard_layer.push_back(epi_layer);
					myocard_layer.push_back(endo_layer);
				
				PorcineCellFactory cell_factory(
					stimulusinfo,
					myocard_layer,									
					meshname,
					&PorcineHeterogen,
					mesh.rGetNodePermutation(),
					ischaemia,
					k0,
					ikatp,
					ina,
					ica,
					gncx,
					pnak,
					inal,
					radius,
					xx,
					yy,
					zz,
					transmural
					); 			
														
				cell_factory.SetHeartGeometryInformation(&heart_geom_info);
			
			
			// Problem creation
				HeartEventHandler::EndEvent(HeartEventHandler::INITIALISE);
				HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
				
				BidomainWithBathProblem<3>* bidomain_problem = new BidomainWithBathProblem<3> ( &cell_factory );			 


				HeartConfig::Instance()->SetOutputDirectory(output_directory);
				HeartConfig::Instance()->SetOutputFilenamePrefix(filename_prefix);
			
				std::vector<unsigned> ECGnodes;
				std::string ecg_nodes_filepath = "/home/user/meshes/" + meshname + "/ecg_nodes.txt";
				std::ifstream ECGnodes_file( ecg_nodes_filepath.c_str() );
				
				while ( ECGnodes_file >> node )
				{
					ECGnodes.push_back(std::atoi(node.c_str()));
				}
				if( !ECGnodes_file.eof() )
				{
					EXCEPTION("Couldn't read ECG_nodes file.");
				}
				
				std::vector<unsigned> groundnodes;
				std::string ground_nodes_filepath = "/home/user/meshes/" + meshname + "/ground_nodes.txt";
				std::ifstream groundnodes_file( ground_nodes_filepath.c_str() );

				while ( groundnodes_file >> node )
				{
					groundnodes.push_back(std::atoi(node.c_str()));
				}
				if( !groundnodes_file.eof() )
				{
					EXCEPTION("Couldn't read ground_nodes file.");
				} 
				
				unsigned j = 0;
				
				boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_LA(new SingleTraceOutputModifier("LA.txt", mesh.rGetNodePermutation()[ECGnodes[j]])); j++;
				bidomain_problem->AddOutputModifier(trace_modifier_LA);
				boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_RA(new SingleTraceOutputModifier("RA.txt", mesh.rGetNodePermutation()[ECGnodes[j]])); j++;
				bidomain_problem->AddOutputModifier(trace_modifier_RA);
				boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_LL(new SingleTraceOutputModifier("LL.txt", mesh.rGetNodePermutation()[ECGnodes[j]])); j++;
				bidomain_problem->AddOutputModifier(trace_modifier_LL);
				boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_RL(new SingleTraceOutputModifier("RL.txt", mesh.rGetNodePermutation()[ECGnodes[j]])); j++;
				bidomain_problem->AddOutputModifier(trace_modifier_RL);
				boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_V1(new SingleTraceOutputModifier("V1.txt", mesh.rGetNodePermutation()[ECGnodes[j]])); j++;
				bidomain_problem->AddOutputModifier(trace_modifier_V1);
				boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_V2a(new SingleTraceOutputModifier("V2.txt", mesh.rGetNodePermutation()[ECGnodes[j]])); j++;
				bidomain_problem->AddOutputModifier(trace_modifier_V2a);
				boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_V3a(new SingleTraceOutputModifier("V3.txt", mesh.rGetNodePermutation()[ECGnodes[j]])); j++;
				bidomain_problem->AddOutputModifier(trace_modifier_V3a);
				boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_V4(new SingleTraceOutputModifier("V4.txt", mesh.rGetNodePermutation()[ECGnodes[j]])); j++;
				bidomain_problem->AddOutputModifier(trace_modifier_V4);
				boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_V5(new SingleTraceOutputModifier("V5.txt", mesh.rGetNodePermutation()[ECGnodes[j]])); j++;
				bidomain_problem->AddOutputModifier(trace_modifier_V5);
				boost::shared_ptr<SingleTraceOutputModifier> trace_modifier_V6(new SingleTraceOutputModifier("V6.txt", mesh.rGetNodePermutation()[ECGnodes[j]]));
				bidomain_problem->AddOutputModifier(trace_modifier_V6);

				if (post_process == true)
				{
					// activation map
			        boost::shared_ptr<ActivationOutputModifier> activation_map(new ActivationOutputModifier("a_map", 0.0));
			        bidomain_problem->AddOutputModifier(activation_map);
			    }

				bidomain_problem->SetMesh( &mesh );
				bidomain_problem->SetWriteInfo();
				bidomain_problem->SetUseHdf5DataWriterCache(true); // test feature

				bidomain_problem->Initialise();

				// PorcineConductivityModifier
					//BidomainWithBathProblem<3>* bidomain_problem = new BidomainWithBathProblem<3> ( &cell_factory );
					//BidomainTissue<3>* p_bidomain_tissue = bidomain_problem->GetBidomainTissue();
					/*
					PorcineConductivityModifier modifier ( 	ischaemia,
															//phase1a,
															//phase1b,
															factor_h,
															radius,
															xx,
															yy,
															zz,
															transmural,
															*bidomain_problem
															//p_bidomain_tissue
														);
					std::cout << "SetConductivityModifier will be set now" << std::endl;
        			//p_bidomain_tissue->SetConductivityModifier( &modifier ); // ok
        			bidomain_problem->GetTissue()->SetConductivityModifier( &modifier );
        			//bidomain_problem->SetConductivityModifier( &modifier );
        			std::cout << "SetConductivityModifier is set now" << std::endl;
        			*/

				std::vector<unsigned> grounded_nodes;											
				grounded_nodes.push_back(mesh.rGetNodePermutation()[groundnodes[j]]); j++;	
				grounded_nodes.push_back(mesh.rGetNodePermutation()[groundnodes[j]]); j++;	
				grounded_nodes.push_back(mesh.rGetNodePermutation()[groundnodes[j]]); j++;
				bidomain_problem->SetFixedExtracellularPotentialNodes(grounded_nodes);

				bidomain_problem->Solve();
			

				HeartEventHandler::Headings();
		        HeartEventHandler::Report();
		        HeartEventHandler::Headings();
		        
		        Warnings::QuietDestroy(); 
		        
				// save permutations
		        OutputFileHandler output_file_handler(output_directory, false);
		        if ( PetscTools::AmMaster() )
		        {
		            out_stream perm_out_stream = output_file_handler.OpenOutputFile("permutations.txt", std::ios::out);
		            std::vector<unsigned> perm_vec = mesh.rGetNodePermutation();
		            for (unsigned i=0; i<perm_vec.size(); ++i)
		            {
		                (*perm_out_stream) << perm_vec[i] << "\n";
		            }
        }
     }
};
