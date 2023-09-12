#include "DifferentialNotchAdhesionForce.hpp"
#include "Cell.hpp"
#include "CollierSrnModel.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
DifferentialNotchAdhesionForce<ELEMENT_DIM, SPACE_DIM>::DifferentialNotchAdhesionForce()
   : GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>(),
     mMaximumSpringConstantMultiplier(1.0),
     mAdhesionDecayConstant(1.0)
{
}


//change this to update the spring multiplier dependents on Notch
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DifferentialNotchAdhesionForce<ELEMENT_DIM, SPACE_DIM>::VariableSpringConstantMultiplicationFactor(
    unsigned nodeAGlobalIndex,
    unsigned nodeBGlobalIndex,
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation,
    bool isCloserThanRestLength)
{
        //Get pointers of cells
        CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
        CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

        //using the pointers, get the Notch values for each of the cells
        CollierSrnModel* p_model_A = static_cast<CollierSrnModel*>(p_cell_A ->GetSrnModel());
        CollierSrnModel* p_model_B = static_cast<CollierSrnModel*>(p_cell_B ->GetSrnModel());

        double Notch_A = p_model_A->GetNotch();
        double Notch_B = p_model_B->GetNotch();

        double delta_Notch_AB = abs(Notch_A - Notch_B);

        //set the spring force multiplier depdent on the notch expression in both cells (decay law)
        return mMaximumSpringConstantMultiplier*exp(-mAdhesionDecayConstant*delta_Notch_AB);

}
//New functions for decay parameters
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DifferentialNotchAdhesionForce<ELEMENT_DIM, SPACE_DIM>::GetMaximumSpringConstantMultiplier()
{
    return mMaximumSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DifferentialNotchAdhesionForce<ELEMENT_DIM, SPACE_DIM>::SetMaximumSpringConstantMultiplier(double typeSpringConstantMultiplier)
{
    assert(typeSpringConstantMultiplier > 0.0);
    mMaximumSpringConstantMultiplier = typeSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DifferentialNotchAdhesionForce<ELEMENT_DIM, SPACE_DIM>::GetAdhesionDecayConstant()
{
    return mAdhesionDecayConstant;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DifferentialNotchAdhesionForce<ELEMENT_DIM, SPACE_DIM>::SetAdhesionDecayConstant(double typeSpringDecay)
{
    assert(typeSpringDecay > 0.0);
    mAdhesionDecayConstant = typeSpringDecay;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DifferentialNotchAdhesionForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MaximumSpringConstantMultiplier>" << mMaximumSpringConstantMultiplier << "</MaximumSpringConstantMultiplier>\n";
    *rParamsFile << "\t\t\t<AdhesionDecayConstants>" << mAdhesionDecayConstant << "</AdhesionDecayConstant>\n";

    // Call direct parent class
    GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class DifferentialNotchAdhesionForce<1,1>;
template class DifferentialNotchAdhesionForce<1,2>;
template class DifferentialNotchAdhesionForce<2,2>;
template class DifferentialNotchAdhesionForce<1,3>;
template class DifferentialNotchAdhesionForce<2,3>;
template class DifferentialNotchAdhesionForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(DifferentialNotchAdhesionForce)