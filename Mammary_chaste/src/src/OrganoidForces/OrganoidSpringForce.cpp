#include "OrganoidSpringForce.hpp"
#include "Cell.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
OrganoidSpringForce<ELEMENT_DIM, SPACE_DIM>::OrganoidSpringForce()
   : GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>(),
     mHomotypicTypeSpringConstantMultiplier(1.0),
     mHeterotypicSpringConstantMultiplier(1.0)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double OrganoidSpringForce<ELEMENT_DIM, SPACE_DIM>::VariableSpringConstantMultiplicationFactor(
    unsigned nodeAGlobalIndex,
    unsigned nodeBGlobalIndex,
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation,
    bool isCloserThanRestLength)
{
    

    if (isCloserThanRestLength)
    {
        return 1.0;
    }
    else
    {
        
        // Determine which (if any) of the cells corresponding to these nodes are labelled
        CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);

        CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);
    
        
        // For heterotypic interactions, scale the spring constant by mHeterotypicSpringConstantMultiplier
        if (p_cell_A->GetCellProliferativeType() != p_cell_B->GetCellProliferativeType() )
        {
            return mHeterotypicSpringConstantMultiplier;
        }
        else 
        {
            // For homotypic interactions between cell of different types, scale the spring constant by mHomotypicLabelledSpringConstantMultiplier
                return mHomotypicTypeSpringConstantMultiplier;
            }
            
            
        }
    }


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double OrganoidSpringForce<ELEMENT_DIM, SPACE_DIM>::GetHomotypicTypeSpringConstantMultiplier()
{
    return mHomotypicTypeSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OrganoidSpringForce<ELEMENT_DIM, SPACE_DIM>::SetHomotypicTypeSpringConstantMultiplier(double typeSpringConstantMultiplier)
{
    assert(typeSpringConstantMultiplier > 0.0);
    mHomotypicTypeSpringConstantMultiplier = typeSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double OrganoidSpringForce<ELEMENT_DIM, SPACE_DIM>::GetHeterotypicSpringConstantMultiplier()
{
    return mHeterotypicSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OrganoidSpringForce<ELEMENT_DIM, SPACE_DIM>::SetHeterotypicSpringConstantMultiplier(double heterotypicSpringConstantMultiplier)
{
    assert(heterotypicSpringConstantMultiplier > 0.0);
    mHeterotypicSpringConstantMultiplier = heterotypicSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OrganoidSpringForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<HomotypicTypeSpringConstantMultiplier>" << mHomotypicTypeSpringConstantMultiplier << "</HomotypicTypeSpringConstantMultiplier>\n";
    *rParamsFile << "\t\t\t<HeterotypicSpringConstantMultiplier>" << mHeterotypicSpringConstantMultiplier << "</HeterotypicSpringConstantMultiplier>\n";

    // Call direct parent class
    GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class OrganoidSpringForce<1,1>;
template class OrganoidSpringForce<1,2>;
template class OrganoidSpringForce<2,2>;
template class OrganoidSpringForce<1,3>;
template class OrganoidSpringForce<2,3>;
template class OrganoidSpringForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OrganoidSpringForce)
