#include "CollierFixedTrackingModifier.hpp"
#include "CollierSrnModel.hpp"
#include "MyoEpiCellProliferativeType.hpp"
#include "BasalStemCellProliferativeType.hpp"
#include "GhostCellProliferativeType.hpp"
#include <AbstractTwoBodyInteractionForce.hpp>
#include "CellLabel.hpp"


template<unsigned DIM>
CollierFixedTrackingModifier<DIM>::CollierFixedTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
    mSameCellW1(0.11)
{
}

template<unsigned DIM>
CollierFixedTrackingModifier<DIM>::~CollierFixedTrackingModifier()
{
}


template<unsigned DIM>
	void CollierFixedTrackingModifier<DIM>::SetSameCellW1(double samecellw1)
	{
	    assert(samecellw1 > 0.0);
	    mSameCellW1 = samecellw1;
	}

template<unsigned DIM>
	double CollierFixedTrackingModifier<DIM>::GetSameCellW1()
	{
	    return mSameCellW1;
	}

template<unsigned DIM>
void CollierFixedTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CollierFixedTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CollierFixedTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();
    double number_of_neighbours = 0;
    int number_of_sweeps = 0;

    //static weight coefficients from connectivity simulations 
    static double w2 = 1;

    //connectivity radius
    static double Connectivity_radius = sqrt(2);
     
    // First recover each cell's Notch and Delta concentrations from the ODEs and store in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        CollierSrnModel* p_model = static_cast<CollierSrnModel*>(cell_iter->GetSrnModel());
        double this_notch = p_model->GetNotch();
        double this_delta = p_model->GetDelta();
      

        // Note that the state variables must be in the same order as listed in MamDeltaNotchOdeSystem
        cell_iter->GetCellData()->SetItem("notch", this_notch);
        cell_iter->GetCellData()->SetItem("delta", this_delta);
       
      
    }
      
    // Next iterate over the population to compute and store each cell's neighbouring Delta concentration in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {   
        //to be overwritten
        std::set<unsigned> neighbour_indices;

       if (bool(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)))
        {   
            assert(static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)->GetMechanicsCutOffLength() > Connectivity_radius);
            unsigned int cell_index = cell_iter->GetCellId(); 
            neighbour_indices = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)->GetNodesWithinNeighbourhoodRadius(cell_index,Connectivity_radius);
        }
        else
        {
            neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);
        }


        // Compute this cell's average neighbouring Delta and Dkk1 concentration and store in CellData
        if (!neighbour_indices.empty())
        {
            //for average neighbour counts
            number_of_neighbours += neighbour_indices.size();
            number_of_sweeps  += 1;
            int num_same_type = 0.0;
            int num_diff_type = 0.0;
            

            double mean_delta = 0.0; //juxtacrine signalling only

            //an initial sweep over neighbours to calculate the cell-type ratio
            for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                 iter != neighbour_indices.end();
                 ++iter)
            {
                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);

                bool neighbour_Cell_is_labelled = p_cell->template HasCellProperty<GhostCellProliferativeType>();
                if (!neighbour_Cell_is_labelled)
                {
                    if (cell_iter->template HasCellProperty<BasalStemCellProliferativeType>() || cell_iter->template HasCellProperty<MyoEpiCellProliferativeType>())
                    {
                        if (p_cell->template HasCellProperty<BasalStemCellProliferativeType>() || p_cell->template HasCellProperty<MyoEpiCellProliferativeType>())
                        {
                            num_same_type += 1;
                        }
                        else
                        {
                            num_diff_type += 1;
                        }
                    }
                    else
                    {
                        if (p_cell->template HasCellProperty<BasalStemCellProliferativeType>() || p_cell->template HasCellProperty<MyoEpiCellProliferativeType>())
                        {
                            num_diff_type += 1;
                        }
                        else
                        {
                            num_same_type += 1;
                        }
                    }
                }
            }


           for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                 iter != neighbour_indices.end();
                 ++iter)
            {

                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);

                bool neighbour_Cell_is_labelled = p_cell->template HasCellProperty<GhostCellProliferativeType>();

                if (!neighbour_Cell_is_labelled)
                {

                    if (cell_iter->template HasCellProperty<BasalStemCellProliferativeType>() || cell_iter->template HasCellProperty<MyoEpiCellProliferativeType>())
                    {
                        if (p_cell->template HasCellProperty<BasalStemCellProliferativeType>() || p_cell->template HasCellProperty<MyoEpiCellProliferativeType>())
                        {
                            double this_delta = mSameCellW1 * (p_cell->GetCellData()->GetItem("delta"));
                            mean_delta += this_delta;
                        }
                        else
                        {
                            double this_delta = w2 * (p_cell->GetCellData()->GetItem("delta"));
                            mean_delta += this_delta;
                        }
                    }
                    else
                    {
                        if (p_cell->template HasCellProperty<BasalStemCellProliferativeType>() || p_cell->template HasCellProperty<MyoEpiCellProliferativeType>())
                        {
                            double this_delta = w2 * (p_cell->GetCellData()->GetItem("delta"));
                            mean_delta += this_delta;
                        }
                        else
                        {
                            double this_delta = mSameCellW1 * (p_cell->GetCellData()->GetItem("delta"));
                            mean_delta += this_delta;
                        }
                    }
                }
            }

             if ((num_same_type*mSameCellW1 + num_diff_type*w2) < 1e-5){
                cell_iter->GetCellData()->SetItem("mean delta", 0);
              }
             else{
                cell_iter->GetCellData()->SetItem("mean delta", mean_delta/(num_same_type*mSameCellW1 + num_diff_type*w2));
               }
        }
        else
        {

            cell_iter->GetCellData()->SetItem("mean delta", 0.0);
        }
        if (cell_iter->template HasCellProperty<GhostCellProliferativeType>()) {
            cell_iter->GetCellData()->SetItem("notch", 0.0);
            cell_iter->GetCellData()->SetItem("delta", 0.0);
        }
       
    }


}

template<unsigned DIM>
void CollierFixedTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class CollierFixedTrackingModifier<1>;
template class CollierFixedTrackingModifier<2>;
template class CollierFixedTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CollierFixedTrackingModifier)