#include "MamDeltaNotchTrackingModifier.hpp"
#include "MamDeltaNotchSrnModel.hpp"
#include "MyoEpiCellProliferativeType.hpp"
#include "BasalStemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include <AbstractTwoBodyInteractionForce.hpp>
#include "CellLabel.hpp"


template<unsigned DIM>
MamDeltaNotchTrackingModifier<DIM>::MamDeltaNotchTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
MamDeltaNotchTrackingModifier<DIM>::~MamDeltaNotchTrackingModifier()
{
}

template<unsigned DIM>
void MamDeltaNotchTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void MamDeltaNotchTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void MamDeltaNotchTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();
    double number_of_neighbours = 0;
    int number_of_sweeps = 0;

    //static weight coefficients from connectivity simulations 
        static double w1 = 0.3;

        static double w2 = 1;

    //connectivity radius
        static double Connectivity_radius = sqrt(2);
     
    // First recover each cell's Notch and Delta concentrations from the ODEs and store in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        MamDeltaNotchSrnModel* p_model = static_cast<MamDeltaNotchSrnModel*>(cell_iter->GetSrnModel());
        double this_notch = p_model->GetNotch();
        double this_hes1 = p_model->GetHes1();
        double this_delta = p_model->GetDelta();
      

        // Note that the state variables must be in the same order as listed in MamDeltaNotchOdeSystem
        cell_iter->GetCellData()->SetItem("notch", this_notch);
        cell_iter->GetCellData()->SetItem("hes1", this_hes1);
        cell_iter->GetCellData()->SetItem("delta", this_delta);
       

       //cout << this_delta << endl;
      
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

            for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                 iter != neighbour_indices.end();
                 ++iter)
            {   

                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                bool neighbour_Cell_is_labelled = p_cell->template HasCellProperty<DifferentiatedCellProliferativeType>();

                //make sure cells in opposing layers are of different types
                if (!neighbour_Cell_is_labelled ){


                    if (p_cell->GetCellProliferativeType() != cell_iter->GetCellProliferativeType() ) {

                        num_diff_type += 1;
                    }

                    else{

                        num_same_type += 1;

                    }

                 }

               
            }


            for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                 iter != neighbour_indices.end();
                 ++iter)
            {   

                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);

                bool neighbour_Cell_is_labelled = p_cell->template HasCellProperty<DifferentiatedCellProliferativeType>();
                if (!neighbour_Cell_is_labelled ){

                    //make sure cells in opposing layers are of different types
                 if (p_cell->GetCellProliferativeType() != cell_iter->GetCellProliferativeType()){

                    double this_delta = w2*(p_cell->GetCellData()->GetItem("delta"));
                    mean_delta += this_delta;
                  
                    }
                    else{ 

                        double this_delta = w1*(p_cell->GetCellData()->GetItem("delta"));
                        mean_delta += this_delta;
      
                    }
                }
 
            }

             if ((num_same_type*w1 + num_diff_type*w2) < 0.0000001){
                cell_iter->GetCellData()->SetItem("mean delta", 0);
              }
             else{
                cell_iter->GetCellData()->SetItem("mean delta", mean_delta/(num_same_type*w1 + num_diff_type*w2));
               }
        }
        else
        {
            // If this cell has no neighbours, such as an isolated cell in a CaBasedCellPopulation, store 0.0 for the cell data
            cell_iter->GetCellData()->SetItem("mean delta", 0.0);
        }

       
    }
    //std::cout<< "Average Neighbour = " <<  number_of_neighbours/number_of_sweeps <<endl;


}

template<unsigned DIM>
void MamDeltaNotchTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class MamDeltaNotchTrackingModifier<1>;
template class MamDeltaNotchTrackingModifier<2>;
template class MamDeltaNotchTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MamDeltaNotchTrackingModifier)