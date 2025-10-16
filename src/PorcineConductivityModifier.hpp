#include "AbstractConductivityModifier.hpp"

//#include "BidomainProblem.hpp" //

class PorcineConductivityModifier : public AbstractConductivityModifier<3,3>
{
private:

    c_matrix<double,3,3> mTensor;

    bool mIschaemia;
    //bool mPhase1a;
    //bool mPhase1b;

    double _factor_h;
	double	_radius;
	double	_xx;
	double	_yy;
	double	_zz;	
	
	unsigned _transmural;

	BidomainProblem<3>& mrProblem; // For tissue, for cells	
	//BidomainTissue<3>& mrProblem;

public:
    PorcineConductivityModifier( 	bool ischaemia,
    								//bool phase1a,
    								//bool phase1b,
    								double factor_h,
									double radius,
									double xx,
									double yy,
									double zz,
									unsigned transmural,
									BidomainProblem<3>& rProblem
									//BidomainTissue<3>& rProblem
									)

        : AbstractConductivityModifier<3,3>(),

        _factor_h(factor_h),
		_radius(radius),
		_xx(xx),
		_yy(yy),
		_zz(zz),
		_transmural(transmural),
		mrProblem(rProblem)

		{
			mIschaemia = ischaemia;
			//mPhase1a = phase1a;
			//mPhase1b = phase1b;
		}

	int counter_inside_matrix =0 ;
	c_matrix<double,3,3>& rCalculateModifiedConductivityTensor(unsigned elementIndex, const c_matrix<double,3,3>& rOriginalConductivity, unsigned domainIndex)		    		
	{
		while (elementIndex < 360000)
		{
		if (counter_inside_matrix == 0)
		{
			std::cout << "got into c_matrix 1st time" << std::endl;
		}
		//int counter_inside_matrix;
		//std::cout << "inside c_matrix number" << counter_inside_matrix << std::endl;
		counter_inside_matrix += 1;

		std::cout << "mIschaemia " << mIschaemia << std::endl;

		if (counter_inside_matrix == 2)
		{
			std::cout << "got into c_matrix 2nd time" << std::endl;
		}

		
		double isch_r = _radius;
        double isch_bz_r = isch_r + 0.5;
        double x_isch = _xx;double y_isch = _yy; double z_isch=_zz;
        ChastePoint<3> isch_centre (x_isch,y_isch,z_isch);
        ChastePoint<3> isch_radius (isch_r,isch_r,isch_r);
        ChastePoint<3> isch_bz_radius (isch_bz_r,isch_bz_r,isch_bz_r);
        ChasteEllipsoid<3> isch_region (isch_centre,isch_radius);
        ChasteEllipsoid<3> isch_bz_region (isch_centre,isch_bz_radius);

        //DistributedVectorFactory* p_factory = mrProblem.rGetMesh().GetDistributedVectorFactory(); // rGetMesh
        //const std::vector<AbstractCardiacCellInterface*>& r_cells = mrProblem.GetTissue()->rGetCellsDistributed();

        //int counter_inside_loop = 0;
        //std::cout << "pre_loop_mofifier" << std::endl;
       
				unsigned i = elementIndex;
				std::cout << "elementIndex" << i << std::endl;
				//std::cout << "counter_inside_loop" << counter_inside_loop << std::endl;
				//std::cout << "inside c_matrix number" << counter_inside_matrix << std::endl;
				//counter_inside_loop += 1;
				
				//
				
				//bool cell_is_in_ischReg = isch_region.DoesContain(mrProblem.GetTissue()->GetNode(i)->rGetLocation()); // rGetMesh
				//GetNode(i) not in AbstractCardiacCellInterface

				//bool cell_is_in_ischReg = isch_region.DoesContain(this->GetTissue()->GetNode(node_index)->rGetLocation()); // GetMesh
		        //bool cell_is_in_ischReg = isch_region.DoesContain(this->rGetMesh()->GetNode(node_index)->rGetLocation()); // GetMesh
		        //bool cell_is_in_isch_BZReg = isch_bz_region.DoesContain(this->GetTissue()->GetNode(node_index)->rGetLocation());

		        //bool cell_is_in_s2StimReg = s2_region.DoesContain(bidomain_problem->rGetMesh().GetNode(i)->rGetLocation());

		        /// x y z
			
		        //double x = this->GetMesh()->GetNode(node_index)->rGetLocation()[0];
		        //double y = this->GetMesh()->GetNode(node_index)->rGetLocation()[1];
		        //double z = this->GetMesh()->GetNode(node_index)->rGetLocation()[2];

		        //double intra_longi_h_n = 1, intra_trans_h_n = 1, intra_normal_h_n = 1;
		        //double intra_longi_h_i = resistance_factor, intra_trans_h_i = resistance_factor, intra_normal_h_i = resistance_factor;

		        double resistance_factor = _factor_h;
		        double intra_h_n = 1;
		        double intra_h_i = resistance_factor;

		        double domain_scaling=1;

		        DistributedVectorFactory* p_factory = mrProblem.rGetMesh().GetDistributedVectorFactory();
			
				for (unsigned k=p_factory->GetLow(); k<p_factory->GetHigh(); k++)
				{
					if (i == k)
					{
						i=k;



		        if (domainIndex == 0) // 0 for intra
					{
						std::cout << "pre_cell_is_in_isch_BZReg?"  << std::endl;
						bool cell_is_in_isch_BZReg = isch_bz_region.DoesContain(mrProblem.rGetMesh().GetNode(i)->rGetLocation()); // rGetMesh
						std::cout << "cell_is_in_isch_BZReg?"  << std::endl;

		        		if (cell_is_in_isch_BZReg && mIschaemia) 
						{
							std::cout << "pre_cell_is_in_ischReg?"  << std::endl;
							bool cell_is_in_ischReg = isch_region.DoesContain(mrProblem.rGetMesh().GetNode(i)->rGetLocation()); // .rGetMesh.
							std::cout << "cell_is_in_ischReg?"  << std::endl;

								if(cell_is_in_ischReg)
								{
									domain_scaling = resistance_factor;	
								}
								else if (_transmural==1)
								{
									double x = mrProblem.rGetMesh().GetNode(i)->rGetLocation()[0]; // rGetMesh
							        double y = mrProblem.rGetMesh().GetNode(i)->rGetLocation()[1];
							        double z = mrProblem.rGetMesh().GetNode(i)->rGetLocation()[2];

									double distance_from_isch_centre = sqrt(pow((x-x_isch),2)+pow((y-y_isch),2)+pow((z-z_isch),2)); // distance from ischaemic centre
									double distance_from_cz = distance_from_isch_centre - isch_r;

									double bz_k = isch_bz_r-isch_r; // K+ border zone width
									double a_intra_h = (intra_h_n-intra_h_i)/bz_k;
								    double y_intra_h = a_intra_h*distance_from_cz+intra_h_i;

									domain_scaling = y_intra_h;
								}
								std::cout << "goingout_cell_is_in_isch_BZReg?_and calc mTensor"  << std::endl;
						}
							//else
					        //{
					        //    domain_scaling = 1; // domainIndex==1 implies extracellular ... assume extra is always constant
					        //}
					}
					else
					{
						domain_scaling = 1;
					}   
					        // save to the "working memory", and return.
							std::cout << "pre_calc_mTensor" << std::endl;
							for ( unsigned j=0; j<1; j++ ) // j=2 no use
							 {
							    mTensor(j,j) = domain_scaling*elementIndex*rOriginalConductivity(j,j);
							 }

							            //mTensor(0,0) = domain_scaling*elementIndex*rOriginalConductivity(0,0); // only intra to modify
							            //std::cout << "mTensor calculated" << std::endl;

							

				}
			}
				
			}
							std::cout << "pre_return_mTensor" << std::endl;
							return mTensor;

					        
	}			
   
};