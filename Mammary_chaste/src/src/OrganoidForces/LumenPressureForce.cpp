#include "LumenPressureForce.hpp"

	template<unsigned DIM>
	LumenPressureForce<DIM>::LumenPressureForce()
	    : AbstractForce<DIM>(),
	          mMovementParameter(0.01),
              mScaleForceParameter(0.3)
	{
	}
	
	template<unsigned DIM>
	LumenPressureForce<DIM>::~LumenPressureForce()
	{
	}
	
template<unsigned DIM>
	void LumenPressureForce<DIM>::SetMovementParameter(double movementParameter)
	{
	    assert(movementParameter > 0.0);
	    mMovementParameter = movementParameter;
	}
	
	template<unsigned DIM>
	double LumenPressureForce<DIM>::GetMovementParameter()
	{
	    return mMovementParameter;
	}


    template<unsigned DIM>
	void LumenPressureForce<DIM>::SetScaleForceParameter(double scaleForceParameter)
	{
	    assert(scaleForceParameter > 0.0);
	    mScaleForceParameter = scaleForceParameter;
	}
	
	template<unsigned DIM>
	double LumenPressureForce<DIM>::GetScaleForceParameter()
	{
	    return mScaleForceParameter;
	}
	
	template<unsigned DIM>
	void LumenPressureForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
	{
	
	
	    // Iterate over the nodes
	    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
	         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
	         ++node_iter)
	    {
	                c_vector<double, DIM> force_contribution;
                    c_vector<double,DIM>  NodeLocation = node_iter->rGetLocation();
                    c_vector<double,DIM> PopCentreLocation = rCellPopulation.GetCentroidOfCellPopulation();
                    double DistFromCentre = norm_2(NodeLocation-PopCentreLocation);
                    c_vector<double,DIM> VectorInOutwardDirection = (1/(DistFromCentre))*(NodeLocation-PopCentreLocation);
	        
            for (unsigned i=0; i<DIM; i++)
	        {
	            /*
	             * force applied in the direction away from the centre of the tissue which is scaled by the distance from the centre
	             */
	
	            force_contribution[i] = mScaleForceParameter*exp(-mMovementParameter*DistFromCentre)*VectorInOutwardDirection[i];
	        }
	        node_iter->AddAppliedForceContribution(force_contribution);
	    }
	}
	
	template<unsigned DIM>
	void LumenPressureForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
	{
	    *rParamsFile << "\t\t\t<MovementParameter>" << mMovementParameter << "</MovementParameter> \n";
	
	    // Call direct parent class
	    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
	}
	
	/////////////////////////////////////////////////////////////////////////////
	// Explicit instantiation
	/////////////////////////////////////////////////////////////////////////////
	
	template class LumenPressureForce<1>;
	template class LumenPressureForce<2>;
	template class LumenPressureForce<3>;
	
	// Serialization for Boost >= 1.36
	#include "SerializationExportWrapperForCpp.hpp"
	EXPORT_TEMPLATE_CLASS_SAME_DIMS(LumenPressureForce)