#include "BscOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

BscOdeSystem::BscOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(8)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<BscOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - Notch concentration for this cell
     * 1 - Hes1 concentration for this cell
     * 2 - Delta concentration for this cell
     * 3 - b-catenin concentration for this cell
     * 4 - TCF/LEF concentration for this cell
     * 5 - TCF Time delay ode
     * 6 - Dickkopf1 concentration for this cell
     * 7 - Neuregulin1 concentration for this cell
     * We store the last state variable so that it can be written
     * to file at each time step alongside the others, and visualized.
     */

    SetDefaultInitialCondition(0, 0.5); // soon overwritten
    SetDefaultInitialCondition(1, 0.5); // soon overwritten
    SetDefaultInitialCondition(2, 0.5); // soon overwritten
    SetDefaultInitialCondition(3, 9.0); // soon overwritten
    SetDefaultInitialCondition(4, 0.5); // soon overwritten
    SetDefaultInitialCondition(5, 0.5); // soon overwritten
    SetDefaultInitialCondition(6, 0.5); // soon overwritten
    SetDefaultInitialCondition(7, 0.5); // soon overwritten

    this->mParameters.push_back(0.5); 
    this->mParameters.push_back(0.5);

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

BscOdeSystem::~BscOdeSystem()
{
}

void BscOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double notch =rY[0];
    double hes1 = rY[1];
    double delta =rY[2];
    double Bcat = rY[3];
    double tcf =  rY[4];
    double tcf_delay = rY[5];
    double dkk1 = rY[6];
    double nrg1 = rY[7];


    double mean_delta = this->mParameters[0]; // Shorthand for "this->mParameter("Mean Delta");"
    double mean_dkk1 = this->mParameters[1];//this->mParameters[1];

    const int r_bar = 1.0; //average rspo1 stimulus from neighbours - later to be implemented from luminal cells
    const int w_bar = 1.0; //average wnt4 stimulus - later to be implemented from luminal cells


    //Definition of the Basal stem cell ODE system.

    rDY[0] = ( 12.8556* pow(mean_delta,3))/(pow(2.0384,3) + pow(mean_delta,3))  - 0.005*notch*Bcat - 1.02*notch; //notch


    rDY[1] = ((1.2*pow(20,3))/(pow(20,3) + pow(hes1,3))) * (  (19.9992* pow(notch,3))/(pow(5.0404,3) + pow(notch,3)) + (0.3555* pow(Bcat,3))/(pow(18.3465,3) + pow(Bcat,3))  ) - 3.9*hes1; //hes1


    rDY[2] = (12.1954*pow(0.6505,3))/(pow(0.6505,3) + pow(hes1,3)) - 2.94*delta; //delta


    rDY[3] = ((30.1695* pow(r_bar,3))/(pow(2.1260,3) + pow(r_bar,3))) * ((30.1695* pow(w_bar,3))/(pow(2.1260,3) + pow(w_bar,3))) * ( (4.2342*pow(1.5561,3))/(pow(1.5561,3) + pow(mean_dkk1,3))) + 3.5604 - 0.005*notch*Bcat - 0.3816*Bcat; //beta-catenin


    rDY[4] = (3.7680* pow(Bcat,3))/(pow(9.2733,3) + pow(Bcat,3)) - 0.3*tcf; //TCF\LEF


    rDY[5] = 0.5*(tcf - tcf_delay); //delay tcf


    rDY[6] = (0.5966* pow(tcf_delay,3))/( pow(5.6760,3) + pow(tcf_delay,3)) - 0.2886*dkk1; //Dickkopf1


    rDY[7] = (33.36* pow(10.0638,3))/( pow(10.0638,3) + pow(notch,3)) - 1.668*nrg1; //Neuregulin1 
}

template<>
void CellwiseOdeSystemInformation<BscOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Notch");
    this->mInitialConditions.push_back(0.5); // will be filled in later

    this->mVariableNames.push_back("Hes1");
    this->mInitialConditions.push_back(0.5); // will be filled in later

    this->mVariableNames.push_back("Delta");
    this->mInitialConditions.push_back(0.5); // will be filled in later

    this->mVariableNames.push_back("Bcat");
    this->mInitialConditions.push_back(9.0); // will be filled in later

    this->mVariableNames.push_back("TCF");
    this->mInitialConditions.push_back(0.5); // will be filled in later

    this->mVariableNames.push_back("TCF_DELAY");
    this->mInitialConditions.push_back(0.5); // will be filled in later

    this->mVariableNames.push_back("Dkk1");
    this->mInitialConditions.push_back(0.5); // will be filled in later

    this->mVariableNames.push_back("Nrg1");
    this->mInitialConditions.push_back(0.5); // will be filled in later

    // If this is ever not the first parameter change the line
    //double mean_delta = this->mParameters[0]; in EvaluateYDerivatives().
    this->mParameterNames.push_back("Mean Delta");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Mean Dkk1");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BscOdeSystem)
