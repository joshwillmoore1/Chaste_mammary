#include "MatureCellKiller.hpp"
#include "CellLabel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "Cell.hpp"
#include "AbstractCellPopulation.hpp"
#include "LumenERPositiveCellProliferativeType.hpp"
#include "LumenERNegativeCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "CellLabel.hpp"


template<unsigned DIM>
MatureCellKiller<DIM>::MatureCellKiller(AbstractCellPopulation<DIM>* pCellPopulation)
: AbstractCellKiller<DIM>(pCellPopulation) 
{
    
}

template<unsigned SPACE_DIM>
void MatureCellKiller<SPACE_DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
    
}

template<unsigned DIM>
void MatureCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath() 
{
    const int age_tol = 80;

     for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {

            if (cell_iter->GetAge() > age_tol)
            {
                cell_iter->StartApoptosis();
            }


    }

}

template<unsigned DIM>
void MatureCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class MatureCellKiller<1>;
template class MatureCellKiller<2>;
template class MatureCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MatureCellKiller)
