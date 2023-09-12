#include "BscSrnModel.hpp"

BscSrnModel::BscSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(8, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<BscSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<BscSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

BscSrnModel::BscSrnModel(const BscSrnModel& rModel)
    : AbstractOdeSrnModel(rModel)
{
    /*
     * Set each member variable of the new SRN model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new SRN model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new SRN model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */

    assert(rModel.GetOdeSystem());
    SetOdeSystem(new BscOdeSystem(rModel.GetOdeSystem()->rGetStateVariables()));
}

AbstractSrnModel* BscSrnModel::CreateSrnModel()
{
    return new BscSrnModel(*this);
}

void BscSrnModel::SimulateToCurrentTime()
{
    // Custom behaviour
    UpdateBsc();

    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void BscSrnModel::Initialise()
{
    AbstractOdeSrnModel::Initialise(new BscOdeSystem);
}

void BscSrnModel::UpdateBsc()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    double mean_delta = mpCell->GetCellData()->GetItem("mean delta");
    mpOdeSystem->SetParameter("Mean Delta", mean_delta);

    double mean_dkk1 = mpCell->GetCellData()->GetItem("mean dkk1");
    mpOdeSystem->SetParameter("Mean Dkk1", mean_dkk1);
}

double BscSrnModel::GetNotch()
{
    assert(mpOdeSystem != nullptr);
    double notch = mpOdeSystem->rGetStateVariables()[0];
    return notch;
}

double BscSrnModel::GetHes1()
{
    assert(mpOdeSystem != nullptr);
    double hes1 = mpOdeSystem->rGetStateVariables()[1];
    return hes1;
}


double BscSrnModel::GetDelta()
{
    assert(mpOdeSystem != nullptr);
    double delta = mpOdeSystem->rGetStateVariables()[2];
    return delta;
}

double BscSrnModel::GetBcat()
{
    assert(mpOdeSystem != nullptr);
    double Bcat = mpOdeSystem->rGetStateVariables()[3];
    return Bcat;
}

double BscSrnModel::GetTCF()
{
    assert(mpOdeSystem != nullptr);
    double tcf = mpOdeSystem->rGetStateVariables()[4];
    return tcf;
}

double BscSrnModel::GetTCF_DELAY()
{
    assert(mpOdeSystem != nullptr);
    double tcf_delay = mpOdeSystem->rGetStateVariables()[5];
    return tcf_delay;
}

double BscSrnModel::GetDkk1()
{
    assert(mpOdeSystem != nullptr);
    double dkk1 = mpOdeSystem->rGetStateVariables()[6];
    return dkk1;
}

double BscSrnModel::GetNrg1()
{
    assert(mpOdeSystem != nullptr);
    double nrg = mpOdeSystem->rGetStateVariables()[7];
    return nrg;
}



double BscSrnModel::GetMeanNeighbouringDelta()
{
    assert(mpOdeSystem != nullptr);
    double mean_neighbouring_delta = mpOdeSystem->GetParameter("Mean Delta");
    return mean_neighbouring_delta;
}

double BscSrnModel::GetMeanNeighbouringDkk1()
{
    assert(mpOdeSystem != nullptr);
    double mean_neighbouring_dkk1 = mpOdeSystem->GetParameter("Mean Dkk1");
    return mean_neighbouring_dkk1;
}

void BscSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BscSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(BscSrnModel)
