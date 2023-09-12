#include "MamDeltaNotchOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

MamDeltaNotchOdeSystem::MamDeltaNotchOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(3)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<MamDeltaNotchOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - Notch concentration for this cell
     * 1 - Hes1 concentration for this cell
     * 2 - Delta concentration for this cell
     * We store the last state variable so that it can be written
     * to file at each time step alongside the others, and visualized.
     */

    SetDefaultInitialCondition(0, 1); // soon overwritten
    SetDefaultInitialCondition(1, 1); // soon overwritten
    SetDefaultInitialCondition(2, 1); // soon overwritten

    this->mParameters.push_back(1);

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

MamDeltaNotchOdeSystem::~MamDeltaNotchOdeSystem()
{
}

void MamDeltaNotchOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double notch =rY[0];
    double hes1 = rY[1];
    double delta =rY[2];
    


    double mean_delta = this->mParameters[0];
     // Shorthand for "this->mParameter("Mean Delta");"

    const double B_cat_ss = 9.0; // we assume here that the crosstalk term is as unstimulated steady state


    //Definition of the NDM ODE system.

    rDY[0] = ( 12.8556* pow(mean_delta,3))/(pow(2.0384,3) + pow(mean_delta,3))  - 0.005*notch*B_cat_ss - 1.02*notch; //notch


    rDY[1] = ((1.2*pow(20,3))/(pow(20,3) + pow(hes1,3))) * (  (19.9992* pow(notch,3))/(pow(5.0404,3) + pow(notch,3)) + (0.3555* pow(B_cat_ss,3))/(pow(18.3465,3) + pow(B_cat_ss,3))  ) - 3.9*hes1; //hes1


    rDY[2] = (12.1954*pow(0.6505,3))/(pow(0.6505,3) + pow(hes1,3)) - 2.94*delta; //delta

}

template<>
void CellwiseOdeSystemInformation<MamDeltaNotchOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Notch");
    this->mInitialConditions.push_back(1); // will be filled in later

    this->mVariableNames.push_back("Hes1");
    this->mInitialConditions.push_back(1); // will be filled in later

    this->mVariableNames.push_back("Delta");
    this->mInitialConditions.push_back(1); // will be filled in later



    // If this is ever not the first parameter change the line
    //double mean_delta = this->mParameters[0]; in EvaluateYDerivatives().
    this->mParameterNames.push_back("Mean Delta");
    this->mParameterUnits.push_back("non-dim");


    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MamDeltaNotchOdeSystem)
