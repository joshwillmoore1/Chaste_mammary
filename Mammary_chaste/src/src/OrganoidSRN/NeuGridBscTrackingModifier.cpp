#include "NeuGridBscTrackingModifier.hpp"
#include "BscSrnModel.hpp"
#include "MyoEpiCellProliferativeType.hpp"
#include "BasalStemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include <AbstractTwoBodyInteractionForce.hpp>


template<unsigned DIM>
NeuGridBscTrackingModifier<DIM>::NeuGridBscTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
NeuGridBscTrackingModifier<DIM>::~NeuGridBscTrackingModifier()
{
}

template<unsigned DIM>
void NeuGridBscTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void NeuGridBscTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void NeuGridBscTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    double number_of_neighbours = 0;
    int number_of_sweeps = 0;

    //static weight coefficients from connectivity simulations
        static double w1 = 0.2;
        static double w2 = 1;

    //connectivity radius
        static double Connectivity_radius = 1;
        static double Length_of_domain = 30;
    


    // First recover each cell's Notch and Delta concentrations from the ODEs and store in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        BscSrnModel* p_model = static_cast<BscSrnModel*>(cell_iter->GetSrnModel());
        double this_notch = p_model->GetNotch();
        double this_hes1 = p_model->GetHes1();
        double this_delta = p_model->GetDelta();
        double this_bcat = p_model->GetBcat();
        double this_tcf = p_model->GetTCF();
        double this_tcf_delay = p_model->GetTCF_DELAY();
        double this_dkk1 = p_model->GetDkk1();
        double this_nrg1 = p_model->GetNrg1();

        // Note that the state variables must be in the same order as listed in BscOdeSystem
        cell_iter->GetCellData()->SetItem("notch", this_notch);
        cell_iter->GetCellData()->SetItem("hes1", this_hes1);
        cell_iter->GetCellData()->SetItem("delta", this_delta);
        cell_iter->GetCellData()->SetItem("Bcat", this_bcat);
        cell_iter->GetCellData()->SetItem("tcf", this_tcf);
        cell_iter->GetCellData()->SetItem("tcf_delay", this_tcf_delay);
        cell_iter->GetCellData()->SetItem("dkk1", this_dkk1);
        cell_iter->GetCellData()->SetItem("nrg1", this_nrg1);
      
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
            /*
            *this needs to be changed if you change the geometry from a neumann grid
            *if the cell is a boundary then place other boundary nodes into set
            */
            if (cell_index == 0){
            neighbour_indices.insert(Length_of_domain -1);
            }
            else if (cell_index == Length_of_domain -1){
            neighbour_indices.insert(0);
            }
            else if (cell_index == Length_of_domain){
            neighbour_indices.insert(2*Length_of_domain-1);
            }
            else if (cell_index == 2*Length_of_domain-1){
            neighbour_indices.insert(Length_of_domain);
            }

        
        }
        else
        {
        std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);
        }
        // Compute this cell's average neighbouring Delta and Dkk1 concentration and store in CellData
        if (!neighbour_indices.empty())
        {
            //for average neighbour counts
            number_of_neighbours += neighbour_indices.size();
            number_of_sweeps  += 1;
            

            double mean_delta = 0.0; //juxtacrine signalling only
            double mean_dkk1 = cell_iter->GetCellData()->GetItem("dkk1"); // autocrine signalling has been demonstrated 

            for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                 iter != neighbour_indices.end();
                 ++iter)
            {   

                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);

                //make sure cells in opposing layers are of different types
                if (p_cell->GetCellProliferativeType() != cell_iter->GetCellProliferativeType()){

                    double this_delta = w2*(p_cell->GetCellData()->GetItem("delta"));
                    mean_delta += this_delta/neighbour_indices.size();

                    double this_dkk1 = p_cell->GetCellData()->GetItem("dkk1");
                    mean_dkk1 += this_dkk1/(neighbour_indices.size()+1);


                }
                else{
                
                    double this_delta = w1*(p_cell->GetCellData()->GetItem("delta"));
                    mean_delta += this_delta/neighbour_indices.size();

                    double this_dkk1 = p_cell->GetCellData()->GetItem("dkk1");
                    mean_dkk1 += this_dkk1/(neighbour_indices.size()+1);
                
                }

               
            }
            cell_iter->GetCellData()->SetItem("mean delta", mean_delta);
            cell_iter->GetCellData()->SetItem("mean dkk1", mean_dkk1);
        }
        
        else
        {
            // If this cell has no neighbours, such as an isolated cell in a CaBasedCellPopulation, store 0.0 for the cell data
            cell_iter->GetCellData()->SetItem("mean delta", 0.0);
            cell_iter->GetCellData()->SetItem("mean dkk1", cell_iter->GetCellData()->GetItem("dkk1")); //autocrine 
        }

       
    }
    std::cout<< "Average Neighbour = " <<  number_of_neighbours/number_of_sweeps <<endl;


}

template<unsigned DIM>
void NeuGridBscTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class NeuGridBscTrackingModifier<1>;
template class NeuGridBscTrackingModifier<2>;
template class NeuGridBscTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NeuGridBscTrackingModifier)