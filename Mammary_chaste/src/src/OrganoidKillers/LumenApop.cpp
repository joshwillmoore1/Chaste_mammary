#include "LumenApop.hpp"
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
LumenApop<DIM>::LumenApop(AbstractCellPopulation<DIM> *pCellPopulation)
    : AbstractCellKiller<DIM>(pCellPopulation)
{
}

template <unsigned SPACE_DIM>
void LumenApop<SPACE_DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
}

template <unsigned DIM>
void LumenApop<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    //const int num_of_lum_neighbours = 4;
    const int JustLumTol = 6;
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

                    if (p_cell->template HasCellProperty<LumenERNegativeCellProliferativeType>() && !p_cell->HasApoptosisBegun() )   {
                        LuminalCellCount++;

                    } else if (p_cell->template HasCellProperty<LumenERPositiveCellProliferativeType>() && !p_cell->HasApoptosisBegun() ){
                        LuminalCellCount++;
                        
                    } 


                    if (p_cell->template HasCellProperty<BasalStemCellProliferativeType>() ||p_cell->template HasCellProperty<MyoEpiCellProliferativeType>() )
                    {
                        isConnectedToBasal = true;
                    }
                }

                //remove the cell if it is surrounded by ghosts
                if (LuminalCellCount >= JustLumTol) 
                {
                    if (cell_iter->GetAge() > AgeTol && isConnectedToBasal == false && !cell_iter->HasApoptosisBegun()){
                        cell_iter->StartApoptosis();
                        this->mpCellPopulation->Update();
                    }
                }

            }

        }
    }

    
}

template <unsigned DIM>
void LumenApop<DIM>::OutputCellKillerParameters(out_stream &rParamsFile)
{

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class LumenApop<1>;
template class LumenApop<2>;
template class LumenApop<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LumenApop)
