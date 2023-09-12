#include "OrganoidG1FixedCellCycleModel.hpp"
#include "BasalStemCellProliferativeType.hpp"
#include "MyoEpiCellProliferativeType.hpp"
#include "GhostCellProliferativeType.hpp"
#include "LumenERPositiveCellProliferativeType.hpp"
#include "LumenERNegativeCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"

OrganoidG1FixedCellCycleModel::OrganoidG1FixedCellCycleModel()
    : AbstractSimplePhaseBasedCellCycleModel()
{
}

OrganoidG1FixedCellCycleModel::OrganoidG1FixedCellCycleModel(const OrganoidG1FixedCellCycleModel &rModel)
    : AbstractSimplePhaseBasedCellCycleModel(rModel)
{
    /*
     * The member variables mGeneration and mMaxTransitGeneration are
     * initialized in the AbstractSimpleGenerationalCellCycleModel
     * constructor.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration is (re)set as soon as InitialiseDaughterCell()
     * is called on the new cell-cycle model.
     * 
     * NOTE: THIS CLASS CURRENTLY HAS ARTIFICIAL CELL CYCLES TIMES -  fix this with Shehata et al. 2018
     * 
     */
}

AbstractCellCycleModel *OrganoidG1FixedCellCycleModel::CreateCellCycleModel()
{
    return new OrganoidG1FixedCellCycleModel(*this);
}

void OrganoidG1FixedCellCycleModel::SetG1Duration()
{

    //Stochastic G1 has been added -  using similar values as the Uniform generational cell cycle model
    RandomNumberGenerator *p_gen = RandomNumberGenerator::Instance();

    //add a reader of generation of the cell

    assert(mpCell != nullptr);

    //organoid cells are set to be diff for SRN investigation

    if (mpCell->GetCellProliferativeType()->IsType<BasalStemCellProliferativeType>())
    {   
        if (p_gen->ranf() > 0.5)
        {
            //mG1Duration = p_gen->GammaRandomDeviate(2, 1.2); this is for the old 
            mG1Duration = p_gen->GammaRandomDeviate(3, 1);
        } else
        {
            mG1Duration = DBL_MAX;
        }
    }
    else if (mpCell->GetCellProliferativeType()->IsType<MyoEpiCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else if (mpCell->GetCellProliferativeType()->IsType<LumenERPositiveCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }

    else if (mpCell->GetCellProliferativeType()->IsType<LumenERNegativeCellProliferativeType>())
    {
        mG1Duration = p_gen->GammaRandomDeviate(5, 1);
    }

    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }

    else if (mpCell->GetCellProliferativeType()->IsType<GhostCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }

    else
    {
        NEVER_REACHED;
    }
}

void OrganoidG1FixedCellCycleModel::SetG2Duration()
{
    //Stochastic G1 has been added -  using similar values as the Uniform generational cell cycle model
    RandomNumberGenerator *p_gen = RandomNumberGenerator::Instance();

    //add a reader of generation of the cell

    assert(mpCell != nullptr);

    //organoid cells are set to be diff for SRN investigation

    if (mpCell->GetCellProliferativeType()->IsType<BasalStemCellProliferativeType>())
    {
        mG2Duration = p_gen->GammaRandomDeviate(3, 1.2);
    }
    else if (mpCell->GetCellProliferativeType()->IsType<MyoEpiCellProliferativeType>())
    {
        mG2Duration = DBL_MAX;
    }
    else if (mpCell->GetCellProliferativeType()->IsType<LumenERPositiveCellProliferativeType>())
    {
        mG2Duration = DBL_MAX;
    }

    else if (mpCell->GetCellProliferativeType()->IsType<LumenERNegativeCellProliferativeType>())
    {
        mG2Duration = p_gen->GammaRandomDeviate(3, 2);
    }

    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG2Duration = DBL_MAX;
    }

    else if (mpCell->GetCellProliferativeType()->IsType<GhostCellProliferativeType>())
    {
        mG2Duration = DBL_MAX;
    }

    else
    {
        NEVER_REACHED;
    }
}

void OrganoidG1FixedCellCycleModel::SetSDuration()
{

    if (mpCell->GetCellProliferativeType()->IsType<GhostCellProliferativeType>())
    {
        mSDuration = DBL_MAX;
    }
    else
    {
        mSDuration = 6;
    }
}

void OrganoidG1FixedCellCycleModel::SetMDuration()
{
    if (mpCell->GetCellProliferativeType()->IsType<GhostCellProliferativeType>())
    {
        mMDuration = DBL_MAX;
    }
    else
    {
        mMDuration = 1;
    }
}

void OrganoidG1FixedCellCycleModel::OutputCellCycleModelParameters(out_stream &rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(OrganoidG1FixedCellCycleModel)