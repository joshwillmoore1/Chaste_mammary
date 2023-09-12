#include "LumenKiller.hpp"
#include "CellLabel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "Cell.hpp"
#include "AbstractCellPopulation.hpp"
#include "LumenERPositiveCellProliferativeType.hpp"
#include "LumenERNegativeCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "GhostCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "CellLabel.hpp"
#include "MyoEpiCellProliferativeType.hpp"
#include "BasalStemCellProliferativeType.hpp"

template <unsigned DIM>
LumenKiller<DIM>::LumenKiller(AbstractCellPopulation<DIM> *pCellPopulation)
    : AbstractCellKiller<DIM>(pCellPopulation)
{
}

template <unsigned SPACE_DIM>
void LumenKiller<SPACE_DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
}

template <unsigned DIM>
void LumenKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    //const int num_of_lum_neighbours = 4;
    const int CombinedGhostsTol = 5;
    const int JustLumTol = 5;
    const double AgeTol = 3;

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
            // only iterate over the live cells
        if ( !(cell_iter->template HasCellProperty<GhostCellProliferativeType>())  )
        {
            // Get the set of neighbouring location indices
            std::set<unsigned> neighbour_indices = this->mpCellPopulation->GetNeighbouringLocationIndices(*cell_iter);
            int GhostCellCount = 0;
            int LuminalCellCount = 0;

            bool isConnectedToBasal = false;

            //if all surrounding cells have labels then kill
            if (!neighbour_indices.empty())
            {

                for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                     iter != neighbour_indices.end();
                     ++iter)
                {
                    CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(*iter);

                    if (p_cell->template HasCellProperty<GhostCellProliferativeType>())
                    {
                      GhostCellCount++;   
                    } else if (p_cell->template HasCellProperty<LumenERNegativeCellProliferativeType>()){
                        LuminalCellCount++;

                    } else if (p_cell->template HasCellProperty<LumenERPositiveCellProliferativeType>()){
                        LuminalCellCount++;
                        
                    } 


                    if (p_cell->template HasCellProperty<BasalStemCellProliferativeType>() ||p_cell->template HasCellProperty<MyoEpiCellProliferativeType>() )
                    {
                        isConnectedToBasal = true;
                    }
                }

                //remove the cell if it is surrounded by ghosts
                if (GhostCellCount == neighbour_indices.size()) 
                {
                    if (cell_iter->GetAge() > AgeTol && isConnectedToBasal == false ){
                        boost::shared_ptr<AbstractCellProperty> p_ghost_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());
                        cell_iter->SetCellProliferativeType(p_ghost_type);
                        cell_iter->GetCellData()->SetItem("target area", 1);
                        cell_iter->SetBirthTime(1e10);

                    }
                } else if (GhostCellCount + LuminalCellCount > CombinedGhostsTol ) 
                {
                    if (cell_iter->GetAge() > AgeTol && isConnectedToBasal == false ){
                        boost::shared_ptr<AbstractCellProperty> p_ghost_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());
                        cell_iter->SetCellProliferativeType(p_ghost_type);
                        cell_iter->GetCellData()->SetItem("target area", 1);
                        cell_iter->SetBirthTime(1e10);

                    }

                } else if (GhostCellCount == 0 && LuminalCellCount > JustLumTol) 
                {
                    if (cell_iter->GetAge() > AgeTol && isConnectedToBasal == false ){
                        boost::shared_ptr<AbstractCellProperty> p_ghost_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());
                        cell_iter->SetCellProliferativeType(p_ghost_type);
                        cell_iter->GetCellData()->SetItem("target area", 1);
                        cell_iter->SetBirthTime(1e10);

                    }
                }

                else if (GhostCellCount >=3 && 2*LuminalCellCount < GhostCellCount) 
                {
                    if (cell_iter->GetAge() > AgeTol && isConnectedToBasal == false ){
                        boost::shared_ptr<AbstractCellProperty> p_ghost_type(CellPropertyRegistry::Instance()->Get<GhostCellProliferativeType>());
                        cell_iter->SetCellProliferativeType(p_ghost_type);
                        cell_iter->GetCellData()->SetItem("target area", 1);
                        cell_iter->SetBirthTime(1e10);

                    }
                }

            }

        }
    }
}

template <unsigned DIM>
void LumenKiller<DIM>::OutputCellKillerParameters(out_stream &rParamsFile)
{

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class LumenKiller<1>;
template class LumenKiller<2>;
template class LumenKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LumenKiller)
