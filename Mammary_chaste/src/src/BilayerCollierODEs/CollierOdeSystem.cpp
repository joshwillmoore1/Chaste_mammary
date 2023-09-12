#include "CollierOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

CollierOdeSystem::CollierOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(2)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<CollierOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - Notch concentration for this cell
     * 1 - Delta concentration for this cell
     * We store the last state variable so that it can be written
     * to file at each time step alongside the others, and visualized.
     */

    SetDefaultInitialCondition(0, 1); // soon overwritten
    SetDefaultInitialCondition(1, 1); // soon overwritten

    this->mParameters.push_back(1);

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

CollierOdeSystem::~CollierOdeSystem()
{
}

void CollierOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double notch =rY[0];
    double delta =rY[1];
    


    double mean_delta = this->mParameters[0];

    //parameter values from Collier et al. 1996

    const double a = 0.01; 
    const double b = 100.0;
    const double mu1 = 1.0;
    const double mu2 = 1.0;
    const double k = 2.0;
    const double h = 1.0;



    //Definition of the Collier lateral-inhibition ODE system.

    rDY[0] = (pow(mean_delta,k))/(a + pow(mean_delta,k)) - mu1*notch;
    rDY[1] = 1/(1+b*pow(notch,h)) - mu2*delta; 
}

template<>
void CellwiseOdeSystemInformation<CollierOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Notch");
    this->mInitialConditions.push_back(1); // will be filled in later

    this->mVariableNames.push_back("Delta");
    this->mInitialConditions.push_back(1); // will be filled in later

    this->mParameterNames.push_back("Mean Delta");
    this->mParameterUnits.push_back("non-dim");


    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CollierOdeSystem)
