#include "BoundaryCellTypeMod.hpp"
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
#include "VertexBasedCellPopulation.hpp"
#include "AbstractCellPopulation.hpp"
#include "MutableVertexMesh.hpp"
#include "VertexElement.hpp"
#include "Node.hpp"

template <unsigned DIM>
BoundaryCellTypeMod<DIM>::BoundaryCellTypeMod(AbstractCellPopulation<DIM> *pCellPopulation)
    : AbstractCellKiller<DIM>(pCellPopulation)
{
}

template <unsigned SPACE_DIM>
void BoundaryCellTypeMod<SPACE_DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
}

template <unsigned DIM>
void BoundaryCellTypeMod<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    //const int num_of_lum_neighbours = 4;
    const int CombinedGhostsTol = 5;
    const int JustLumTol = 5;
    const double AgeTol = 3;

    //assert that this can only be used for vertex simulations
    assert(dynamic_cast<VertexBasedCellPopulation<DIM> *>(this->mpCellPopulation));
    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM> *p_cell_population = static_cast<VertexBasedCellPopulation<DIM> *>(this->mpCellPopulation);
    
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {

            /*IMMEDIATE SET OF CELL TYPE BASED ON POSITION IN THE TISSUE
                * if the cell is a boundary cell then it is basal
                * else it is luminal
                * then allow for the ghost cell assignment
                * 
                */

               //get the nodes associated with the cell
                VertexElement<DIM, DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter);
                 if (p_element->IsElementOnBoundary() && cell_iter->GetAge() > AgeTol)
                 {
                     boost::shared_ptr<AbstractCellProperty> p_basal_type(CellPropertyRegistry::Instance()->Get<BasalStemCellProliferativeType>());
                     cell_iter->SetCellProliferativeType(p_basal_type);
                        cell_iter->GetCellData()->SetItem("target area", 0.85);
                 }




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


                    if (p_cell->template HasCellProperty<BasalStemCellProliferativeType>() || p_cell->template HasCellProperty<MyoEpiCellProliferativeType>() )
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
void BoundaryCellTypeMod<DIM>::OutputCellKillerParameters(out_stream &rParamsFile)
{

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class BoundaryCellTypeMod<1>;
template class BoundaryCellTypeMod<2>;
template class BoundaryCellTypeMod<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BoundaryCellTypeMod)