#include "BoundaryCurvatureCalculation.hpp"
#include <AbstractTwoBodyInteractionForce.hpp>
#include "Cell.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "AbstractCellPopulation.hpp"
#include "MutableVertexMesh.hpp"
#include "VertexElement.hpp"
#include "Node.hpp"
#include "SimulationTime.hpp"

template <unsigned DIM>
BoundaryCurvatureCalculation<DIM>::BoundaryCurvatureCalculation()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mSimulationEndTime(100),
      mSimulationDt(0.001),
      mCircleRadius(10),
      mOscRadius(0),
      mBoundarySteps(10)
{
    c_vector<double, DIM> tempCentre;
    for (int i = 0; i < tempCentre.size(); i++)
    {
        tempCentre[i] = 0;
    }
    mCentre = tempCentre;
}

template <unsigned DIM>
BoundaryCurvatureCalculation<DIM>::~BoundaryCurvatureCalculation()
{
}

template <unsigned DIM>
void BoundaryCurvatureCalculation<DIM>::SetSimulationEndTime(double simulationEndTime)
{
    assert(simulationEndTime > 0.0);
    mSimulationEndTime = simulationEndTime;
}

template <unsigned DIM>
double BoundaryCurvatureCalculation<DIM>::GetSimulationEndTime()
{
    return mSimulationEndTime;
}

template <unsigned DIM>
void BoundaryCurvatureCalculation<DIM>::SetSimulationDt(double simulationDt)
{
    assert(simulationDt > 0.0);
    mSimulationDt = simulationDt;
}

template <unsigned DIM>
double BoundaryCurvatureCalculation<DIM>::GetSimulationDt()
{
    return mSimulationDt;
}

template <unsigned DIM>
void BoundaryCurvatureCalculation<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM> &rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

//New get and set functions for the curvature calculations

template <unsigned DIM>
double BoundaryCurvatureCalculation<DIM>::GetCircleRadius()
{
    return mCircleRadius;
}

template <unsigned DIM>
void BoundaryCurvatureCalculation<DIM>::SetCircleRadius(double circleRad)
{
    assert(circleRad > 0.0);
    mCircleRadius = circleRad;
}

template <unsigned DIM>
double BoundaryCurvatureCalculation<DIM>::GetOscRadius()
{
    return mOscRadius;
}

template <unsigned DIM>
void BoundaryCurvatureCalculation<DIM>::SetOscRadius(double oscRad)
{
    assert(oscRad > 0.0);
    mOscRadius = oscRad;
}

template <unsigned DIM>
c_vector<double, DIM> BoundaryCurvatureCalculation<DIM>::GetCentre()
{
    return mCentre;
}

template <unsigned DIM>
void BoundaryCurvatureCalculation<DIM>::SetCentre(c_vector<double, DIM> centre)
{
    mCentre = centre;
}

template <unsigned DIM>
int BoundaryCurvatureCalculation<DIM>::GetBoundarySteps()
{
    return mBoundarySteps;
}

template <unsigned DIM>
void BoundaryCurvatureCalculation<DIM>::SetBoundarySteps(int boundSteps)
{
    assert(boundSteps > 0.0);
    mBoundarySteps = boundSteps;
}

template <unsigned DIM>
void BoundaryCurvatureCalculation<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM> &rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template <unsigned DIM>
void BoundaryCurvatureCalculation<DIM>::UpdateCellData(AbstractCellPopulation<DIM, DIM> &rCellPopulation)
{

    //this should only occur at the final step
    double currentSimTim = SimulationTime::Instance()->GetTime();
    if (currentSimTim == (mSimulationEndTime))
    {
        std::vector<double> cellPositionX;
        std::vector<double> cellPositionY;

        static double DistScale = 7.92;

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            cellPositionX.push_back(cell_location[0]);
            cellPositionY.push_back(cell_location[1]);
        }

        const int arcSteps = 1000;
        const int nOsc = 4;
        //boundaryLinspace
        c_vector<double, arcSteps + 1> BoundaryPointX0;
        c_vector<double, arcSteps + 1> BoundaryPointY0;

        c_vector<double, arcSteps + 1> Curvature;
        c_vector<int, arcSteps + 1> ClosestCellId;

        std::vector<double> BoundaryTheta(arcSteps);
        //c_vector<double, 3*arcSteps> BoundaryThetaLong =  BoundaryThetaLong.
        c_vector<double, arcSteps> BoundaryCurvatures;

        for (double i = 0; i < BoundaryTheta.size(); i++)
        {
            int Index = (int)round(i);
            BoundaryTheta[Index] = 2 * 3.14159265359 * (i / 1000);
        }
        //make a buffer of theta for plus and minus
        std::vector<double> BoundaryThetaLong;
        BoundaryThetaLong.insert(BoundaryThetaLong.end(), BoundaryTheta.begin(), BoundaryTheta.end());
        BoundaryThetaLong.insert(BoundaryThetaLong.end(), BoundaryTheta.begin(), BoundaryTheta.end());
        BoundaryThetaLong.insert(BoundaryThetaLong.end(), BoundaryTheta.begin(), BoundaryTheta.end());

        for (int i = 0; i < (arcSteps + 1); i++)
        {
            double Theta0 = BoundaryThetaLong[arcSteps + i];
            double Theta1 = BoundaryThetaLong[arcSteps + i - mBoundarySteps];
            double Theta2 = BoundaryThetaLong[arcSteps + i + mBoundarySteps];

            double Rad0 = mCircleRadius + mOscRadius * sin(nOsc * Theta0);
            double Rad1 = mCircleRadius + mOscRadius * sin(nOsc * Theta1);
            double Rad2 = mCircleRadius + mOscRadius * sin(nOsc * Theta2);

            double Rad1p = mOscRadius * nOsc * cos(nOsc * Theta1);
            double Rad2p = mOscRadius * nOsc * cos(nOsc * Theta2);

            //for reference
            BoundaryPointX0[i] = Rad0 * cos(Theta0) + mCentre[0];
            BoundaryPointY0[i] = Rad0 * sin(Theta0) + mCentre[1];

            // boundary points either side

            //POINT 1
            double BoundaryPointX1 = Rad1 * cos(Theta1);
            double BoundaryPointY1 = Rad1 * sin(Theta1);

            double BoundaryPointX1p = Rad1p * cos(Theta1) - Rad1 * sin(Theta1);
            double BoundaryPointY1p = Rad1p * sin(Theta1) + Rad1 * cos(Theta1);

            c_vector<double, 2> TangentNormVec1;
            TangentNormVec1[0] = BoundaryPointX1p;
            TangentNormVec1[1] = BoundaryPointY1p;

            double BoundaryPointX1pN = BoundaryPointX1p / norm_2(TangentNormVec1);
            double BoundaryPointY1pN = BoundaryPointY1p / norm_2(TangentNormVec1);

            double Grad1 = BoundaryPointY1pN / BoundaryPointX1pN;

            //POINT 2
            double BoundaryPointX2 = Rad2 * cos(Theta2);
            double BoundaryPointY2 = Rad2 * sin(Theta2);

            double BoundaryPointX2p = Rad2p * cos(Theta2) - Rad2 * sin(Theta2);
            double BoundaryPointY2p = Rad2p * sin(Theta2) + Rad2 * cos(Theta2);

            c_vector<double, 2> TangentNormVec2;
            TangentNormVec2[0] = BoundaryPointX2p;
            TangentNormVec2[1] = BoundaryPointY2p;

            double BoundaryPointX2pN = BoundaryPointX2p / norm_2(TangentNormVec2);
            double BoundaryPointY2pN = BoundaryPointY2p / norm_2(TangentNormVec2);

            double Grad2 = BoundaryPointY2pN / BoundaryPointX2pN;

            //Make sure tangents aren't parallel
            if (abs(Grad1 - Grad2) > 1e-5)
            {
                //get normals and linear coefficients
                double NormalX1 = -BoundaryPointY1pN + BoundaryPointX1;
                double NormalY1 = BoundaryPointX1pN + BoundaryPointY1;

                double M1 = (NormalY1 - BoundaryPointY1) / (NormalX1 - BoundaryPointX1);
                double C1 = NormalY1 - M1 * NormalX1;

                double NormalX2 = -BoundaryPointY2pN + BoundaryPointX2;
                double NormalY2 = BoundaryPointX2pN + BoundaryPointY2;

                double M2 = (NormalY2 - BoundaryPointY2) / (NormalX2 - BoundaryPointX2);
                double C2 = NormalY2 - M2 * NormalX2;

                //get the intersection of the linear equations
                double CentreOfCirX = (C2 - C1) / (M1 - M2);
                double CentreOfCirY = M1 * CentreOfCirX + C1;

                c_vector<double, 2> RadiusOfCirVec;
                RadiusOfCirVec[0] = CentreOfCirX - BoundaryPointX1;
                RadiusOfCirVec[1] = CentreOfCirY - BoundaryPointY1;

                double RadiusOfCirc = norm_2(RadiusOfCirVec) * DistScale;

                //determine whether the curvature is positive or negative

                c_vector<double, 2> MidpointVec;
                MidpointVec[0] = (BoundaryPointX1 + BoundaryPointX2) / 2;
                MidpointVec[1] = (BoundaryPointY1 + BoundaryPointY2) / 2;

                double AngleOfMidpoint = atan2(MidpointVec[1], MidpointVec[0]);
                double TempRad = mCircleRadius + mOscRadius * sin(nOsc * AngleOfMidpoint);

                //possible error here as we haven't shifted for the centre?
                if (TempRad > norm_2(MidpointVec))
                {
                    Curvature[i] = 1 / RadiusOfCirc;
                }
                else
                {
                    Curvature[i] = -1 / RadiusOfCirc;
                }
            }
            else
            {
                Curvature[i] = 0;
            }

            //std::cout << Curvature[i] << std::endl;
            double tempDist = 1e5;
            for (int k = 0; k < cellPositionX.size(); k++)
            {
                c_vector<double, 2> BoundaryPointVec;
                BoundaryPointVec[0] = BoundaryPointX0[i];
                BoundaryPointVec[1] = BoundaryPointY0[i];

                c_vector<double, 2> CellVec;
                CellVec[0] = cellPositionX[k];
                CellVec[1] = cellPositionY[k];

                double MeasureDist = norm_2(BoundaryPointVec - CellVec);
                if (MeasureDist < tempDist)
                {
                    tempDist = MeasureDist;
                    ClosestCellId[i] = k;
                }
            }

            //now find the closest cell and add this curvature to the list then take an average
        }
        int cellIterator = 0;
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            std::vector<double> thisCurvature;

            for (int i = 0; i < (arcSteps + 1); i++)
            {
                if (ClosestCellId[i] == cellIterator)
                {
                    thisCurvature.push_back(Curvature[i]);
                }
            }
            // if the cell is not close to the boundary then it is assigned a weird number to be ignored
            if (thisCurvature.empty())
            {
                cell_iter->GetCellData()->SetItem("Curvature", -100);
            }
            else
            {
                sort(thisCurvature.begin(), thisCurvature.end());
                int size_of_vector = thisCurvature.size();
                double MedianVal = -100;
                if (size_of_vector % 2 == 0)
                {
                    int medIndex = size_of_vector / 2;
                    MedianVal = thisCurvature[medIndex];
                }
                else
                {
                    int medIndexM = floor(size_of_vector * 0.5);
                    int medIndexP = ceil(size_of_vector * 0.5);
                    MedianVal = (thisCurvature[medIndexM] + thisCurvature[medIndexP]) * 0.5;
                }

                cell_iter->GetCellData()->SetItem("Curvature", MedianVal);
            }
            cellIterator++;
        }

        //finally iterate through the cells and count the connected cells with and without curvatures
       
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            double ThisCellCurvature = cell_iter->GetCellData()->GetItem("Curvature");

            if (ThisCellCurvature == -100)
            {
                cell_iter->GetCellData()->SetItem("BoundaryNeigh", -10);
                cell_iter->GetCellData()->SetItem("NonBoundaryNeigh", -10);
            }
            else
            {
                //get neighbours

                std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);
                int numBound = 0;
                int numNonBound = 0;
                if (!neighbour_indices.empty())
                {
                    for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                         iter != neighbour_indices.end();
                         ++iter)
                    {

                            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                            double ThatCellCurvature = p_cell->GetCellData()->GetItem("Curvature");

                            if (ThatCellCurvature == -100)
                            {
                                numNonBound++;
                            } else 
                            {
                                numBound++;
                            }
                    }
                }

                cell_iter->GetCellData()->SetItem("BoundaryNeigh", numBound);
                cell_iter->GetCellData()->SetItem("NonBoundaryNeigh", numNonBound);

            }
        }
    }

    rCellPopulation.Update();
}

template <unsigned DIM>
void BoundaryCurvatureCalculation<DIM>::OutputSimulationModifierParameters(out_stream &rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class BoundaryCurvatureCalculation<1>;
template class BoundaryCurvatureCalculation<2>;
template class BoundaryCurvatureCalculation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BoundaryCurvatureCalculation)