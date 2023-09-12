#include "NagaiHondaDiffAdTypeDependent.hpp"
#include "LumenERPositiveCellProliferativeType.hpp"
#include "LumenERNegativeCellProliferativeType.hpp"
#include "BasalStemCellProliferativeType.hpp"
#include "GhostCellProliferativeType.hpp"
#include "MyoEpiCellProliferativeType.hpp"

template <unsigned DIM>
NagaiHondaDiffAdTypeDependent<DIM>::NagaiHondaDiffAdTypeDependent()
    : NagaiHondaForce<DIM>(),
      mNagaiHondaBasalLuminalAdhesionEnergyParameter(1.0),
      mNagaiHondaLuminalLuminalAdhesionEnergyParameter(1.0),
      mNagaiHondaBasalBasalAdhesionEnergyParameter(1.0),
      mNagaiHondaGhostGhostAdhesionEnergyParameter(1.0),
      mNagaiHondaBasalBoundaryAdhesionEnergyParameter(1.0),
      mNagaiHondaLuminalBoundaryAdhesionEnergyParameter(1.0),
      mNagaiHondaGhostBoundaryAdhesionEnergyParameter(1.0),
      mNagaiHondaLuminalGhostAdhesionEnergyParameter(1.0),
      mNagaiHondaBasalGhostAdhesionEnergyParameter(1.0)
{
}

template <unsigned DIM>
NagaiHondaDiffAdTypeDependent<DIM>::~NagaiHondaDiffAdTypeDependent()
{
}

template <unsigned DIM>
double NagaiHondaDiffAdTypeDependent<DIM>::GetAdhesionParameter(Node<DIM> *pNodeA,
                                                                Node<DIM> *pNodeB,
                                                                VertexBasedCellPopulation<DIM> &rVertexCellPopulation)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        unsigned element_index = *(shared_elements.begin());

        // Get cell associated with this element
        CellPtr p_cell = rVertexCellPopulation.GetCellUsingLocationIndex(element_index);

        // BOUNDARY CELLS ADHESION

        if (p_cell->template HasCellProperty<BasalStemCellProliferativeType>() || p_cell->template HasCellProperty<MyoEpiCellProliferativeType>())
        {
            // This cell is labelled
            return this->GetNagaiHondaBasalBoundaryAdhesionEnergyParameter();
        }
        else if (p_cell->template HasCellProperty<LumenERPositiveCellProliferativeType>() || p_cell->template HasCellProperty<LumenERNegativeCellProliferativeType>())
        {
            // This cell is not labelled
            return this->GetNagaiHondaLuminalBoundaryAdhesionEnergyParameter();
        }
        else if (p_cell->template HasCellProperty<GhostCellProliferativeType>())
        {
            // This cell is not labelled
            return this->GetNagaiHondaGhostBoundaryAdhesionEnergyParameter();
        }
        else
        {
            NEVER_REACHED;
        }
    }
    else
    {
        // Work out the number of labelled cells: 0,1 or 2

        //HERE IS WHERE WE CHANGE FOR DIFFERENT CELL CONNECTIONS
        unsigned num_basal_cells = 0;
        unsigned num_lum_cells = 0;
        unsigned num_ghost_cells = 0;
        for (std::set<unsigned>::iterator iter = shared_elements.begin();
             iter != shared_elements.end();
             ++iter)
        {
            unsigned element_index = *(iter);

            // Get cell associated with this element
            CellPtr p_cell = rVertexCellPopulation.GetCellUsingLocationIndex(element_index);

            if (p_cell->template HasCellProperty<BasalStemCellProliferativeType>() || p_cell->template HasCellProperty<MyoEpiCellProliferativeType>())
            {
                num_basal_cells++;
            }

            else if (p_cell->template HasCellProperty<LumenERPositiveCellProliferativeType>() || p_cell->template HasCellProperty<LumenERNegativeCellProliferativeType>())
            {
                num_lum_cells++;
            }
            else if (p_cell->template HasCellProperty<GhostCellProliferativeType>())
            {
                num_ghost_cells++;
            }
            else
            {
                NEVER_REACHED;
            }
        }

        if (num_basal_cells == 2)
        {
            // Both cells are labelled
            return this->GetNagaiHondaBasalBasalAdhesionEnergyParameter();
        }
        else if (num_lum_cells == 2)
        {
            // One cell is labelled
            return this->GetNagaiHondaLuminalLuminalAdhesionEnergyParameter();
        }
        else if (num_ghost_cells == 2)
        {
            // Neither cell is labelled
            return this->GetNagaiHondaGhostGhostAdhesionEnergyParameter();
        }

        else if (num_ghost_cells == 1 && num_basal_cells == 1)
        {
            // Neither cell is labelled
            return this->GetNagaiHondaBasalGhostAdhesionEnergyParameter();
        }

        else if (num_ghost_cells == 1 && num_lum_cells == 1)
        {
            // Neither cell is labelled
            return this->GetNagaiHondaLuminalGhostAdhesionEnergyParameter();
        }

        else if (num_basal_cells == 1 && num_lum_cells == 1)
        {
            // Neither cell is labelled
            return this->GetNagaiHondaBasalLuminalAdhesionEnergyParameter();
        }
        else
        {
            NEVER_REACHED;
        }
    }
}


/*
*Here is where I am adding the new parameter get functions
*/

template <unsigned DIM>
double NagaiHondaDiffAdTypeDependent<DIM>::GetNagaiHondaBasalLuminalAdhesionEnergyParameter()
{
    return mNagaiHondaBasalLuminalAdhesionEnergyParameter;
}

template <unsigned DIM>
double NagaiHondaDiffAdTypeDependent<DIM>::GetNagaiHondaLuminalLuminalAdhesionEnergyParameter()
{
    return mNagaiHondaLuminalLuminalAdhesionEnergyParameter;
}

template <unsigned DIM>
double NagaiHondaDiffAdTypeDependent<DIM>::GetNagaiHondaBasalBasalAdhesionEnergyParameter()
{
    return mNagaiHondaBasalBasalAdhesionEnergyParameter;
}

template <unsigned DIM>
double NagaiHondaDiffAdTypeDependent<DIM>::GetNagaiHondaGhostGhostAdhesionEnergyParameter()
{
    return mNagaiHondaGhostGhostAdhesionEnergyParameter;
}

template <unsigned DIM>
double NagaiHondaDiffAdTypeDependent<DIM>::GetNagaiHondaBasalBoundaryAdhesionEnergyParameter()
{
    return mNagaiHondaBasalBoundaryAdhesionEnergyParameter;
}

template <unsigned DIM>
double NagaiHondaDiffAdTypeDependent<DIM>::GetNagaiHondaLuminalBoundaryAdhesionEnergyParameter()
{
    return mNagaiHondaLuminalBoundaryAdhesionEnergyParameter;
}

template <unsigned DIM>
double NagaiHondaDiffAdTypeDependent<DIM>::GetNagaiHondaGhostBoundaryAdhesionEnergyParameter()
{
    return mNagaiHondaGhostBoundaryAdhesionEnergyParameter;
}

template <unsigned DIM>
double NagaiHondaDiffAdTypeDependent<DIM>::GetNagaiHondaLuminalGhostAdhesionEnergyParameter()
{
    return mNagaiHondaLuminalGhostAdhesionEnergyParameter;
}

template <unsigned DIM>
double NagaiHondaDiffAdTypeDependent<DIM>::GetNagaiHondaBasalGhostAdhesionEnergyParameter()
{
    return mNagaiHondaBasalGhostAdhesionEnergyParameter;
}

/*
now for the set functions
*/

template <unsigned DIM>
void NagaiHondaDiffAdTypeDependent<DIM>::SetNagaiHondaBasalLuminalAdhesionEnergyParameter(double NagaiHondaBasalLuminalAdhesionEnergyParameter)
{
    mNagaiHondaBasalLuminalAdhesionEnergyParameter = NagaiHondaBasalLuminalAdhesionEnergyParameter;
}

template <unsigned DIM>
void NagaiHondaDiffAdTypeDependent<DIM>::SetNagaiHondaLuminalLuminalAdhesionEnergyParameter(double NagaiHondaLuminalLuminalAdhesionEnergyParameter)
{
    mNagaiHondaLuminalLuminalAdhesionEnergyParameter = NagaiHondaLuminalLuminalAdhesionEnergyParameter;
}

template <unsigned DIM>
void NagaiHondaDiffAdTypeDependent<DIM>::SetNagaiHondaBasalBasalAdhesionEnergyParameter(double NagaiHondaBasalBasalAdhesionEnergyParameter)
{
    mNagaiHondaBasalBasalAdhesionEnergyParameter = NagaiHondaBasalBasalAdhesionEnergyParameter;
}

template <unsigned DIM>
void NagaiHondaDiffAdTypeDependent<DIM>::SetNagaiHondaGhostGhostAdhesionEnergyParameter(double NagaiHondaGhostGhostAdhesionEnergyParameter)
{
    mNagaiHondaGhostGhostAdhesionEnergyParameter = NagaiHondaGhostGhostAdhesionEnergyParameter;
}

template <unsigned DIM>
void NagaiHondaDiffAdTypeDependent<DIM>::SetNagaiHondaBasalBoundaryAdhesionEnergyParameter(double NagaiHondaBasalBoundaryAdhesionEnergyParameter)
{
    mNagaiHondaBasalBoundaryAdhesionEnergyParameter = NagaiHondaBasalBoundaryAdhesionEnergyParameter;
}

template <unsigned DIM>
void NagaiHondaDiffAdTypeDependent<DIM>::SetNagaiHondaLuminalBoundaryAdhesionEnergyParameter(double NagaiHondaLuminalBoundaryAdhesionEnergyParameter)
{
    mNagaiHondaLuminalBoundaryAdhesionEnergyParameter = NagaiHondaLuminalBoundaryAdhesionEnergyParameter;
}

template <unsigned DIM>
void NagaiHondaDiffAdTypeDependent<DIM>::SetNagaiHondaGhostBoundaryAdhesionEnergyParameter(double NagaiHondaGhostBoundaryAdhesionEnergyParameter)
{
    mNagaiHondaGhostBoundaryAdhesionEnergyParameter = NagaiHondaGhostBoundaryAdhesionEnergyParameter;
}

template <unsigned DIM>
void NagaiHondaDiffAdTypeDependent<DIM>::SetNagaiHondaLuminalGhostAdhesionEnergyParameter(double NagaiHondaLuminalGhostAdhesionEnergyParameter)
{
    mNagaiHondaLuminalGhostAdhesionEnergyParameter = NagaiHondaLuminalGhostAdhesionEnergyParameter;
}

template <unsigned DIM>
void NagaiHondaDiffAdTypeDependent<DIM>::SetNagaiHondaBasalGhostAdhesionEnergyParameter(double NagaiHondaBasalGhostAdhesionEnergyParameter)
{
    mNagaiHondaBasalGhostAdhesionEnergyParameter = NagaiHondaBasalGhostAdhesionEnergyParameter;
}

// Explicit instantiation
template class NagaiHondaDiffAdTypeDependent<1>;
template class NagaiHondaDiffAdTypeDependent<2>;
template class NagaiHondaDiffAdTypeDependent<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NagaiHondaDiffAdTypeDependent)
