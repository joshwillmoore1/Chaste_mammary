/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "VertexBilayerBoundaryForce.hpp"
#include "MathsCustomFunctions.hpp"

template <unsigned DIM>
VertexBilayerBoundaryForce<DIM>::VertexBilayerBoundaryForce(double forceStrength, double circleRadius, double oscRadius, c_vector<double, DIM> CentreOfTissue)
    : AbstractForce<DIM>(),
      mForceStrength(forceStrength),
      mCircleRadius(circleRadius),
      mOscRadius(oscRadius),
      mCentreOfTissue(CentreOfTissue)
{
    // We don't want the force to act in the wrong direction
    assert(mForceStrength > 0.0);
    assert(mCircleRadius > 0.0);
    assert(mOscRadius > 0.0);
}

template <unsigned DIM>
VertexBilayerBoundaryForce<DIM>::~VertexBilayerBoundaryForce()
{
}

template <unsigned DIM>
void VertexBilayerBoundaryForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM> &rCellPopulation)
{
    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM> *p_cell_population = static_cast<VertexBasedCellPopulation<DIM> *>(&rCellPopulation);

    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM> *>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("VertexBilayerBoundaryForce is to be used with VertexBasedCellPopulations only");
    }

    c_vector<double, DIM> ReferenceVector;
    ReferenceVector[0] = 0;
    ReferenceVector[1] = 1;

    //boundaryLinspace
    c_vector<double, 1000> BoundaryPointX;
    c_vector<double, 1000> BoundaryPointY;
    for (double i = 0; i < 1000; i++)
    {
        double tempT = 2 * 3.14159265359 * (i / 1000);
        double TempR = mCircleRadius + mOscRadius * sin(4 * tempT);
        int Index = (int)round(i);

        BoundaryPointX[Index] = TempR * cos(tempT) + mCentreOfTissue[0];
        BoundaryPointY[Index] = TempR * sin(tempT) + mCentreOfTissue[1];
    }

    // Iterate over nodes
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = p_cell_population->rGetMesh().GetNodeIteratorBegin();
         node_iter != p_cell_population->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        c_vector<double, DIM> boundary_force = zero_vector<double>(DIM);
        c_vector<double, DIM> NodeLocation = node_iter->rGetLocation();
        c_vector<double, DIM> NodeLocation_trans = NodeLocation - mCentreOfTissue;

        double DistToCentre = norm_2(NodeLocation_trans);

        if (DistToCentre > 1e-4)
        {
            double theta = atan2(NodeLocation_trans[1],NodeLocation_trans[0]);
            double OscilatingRad = mCircleRadius + mOscRadius * sin(4 * theta);

            if (DistToCentre > OscilatingRad)
            {

                c_vector<double, DIM> BoundaryPoint;
                c_vector<double, DIM> nearest_point;
                double SmallestDist = 1e10;
                
                for (int i = 0; i < 1000; i++)
                {

                    BoundaryPoint[0] = BoundaryPointX[i];
                    BoundaryPoint[1] = BoundaryPointY[i];

                    double TempDist = norm_2(BoundaryPoint - NodeLocation);

                    if (TempDist < SmallestDist)
                    {
                        SmallestDist = TempDist;
                        nearest_point = BoundaryPoint;
                    }

                }
                c_vector<double, DIM> InwardVector = nearest_point - NodeLocation;

                for (int i = 0; i < BoundaryPoint.size(); ++i)
                {
                    boundary_force[i] = InwardVector[i] * (mForceStrength)*pow(norm_2(InwardVector), 2);
                }

                node_iter->AddAppliedForceContribution(boundary_force);
            }
        }
    }
}

template <unsigned DIM>
double VertexBilayerBoundaryForce<DIM>::GetForceStrength() const
{
    return mForceStrength;
}

template <unsigned DIM>
void VertexBilayerBoundaryForce<DIM>::OutputForceParameters(out_stream &rParamsFile)
{
    *rParamsFile << "\t\t\t<ForceStrength>" << mForceStrength << "</ForceStrength>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class VertexBilayerBoundaryForce<1>;
template class VertexBilayerBoundaryForce<2>;
template class VertexBilayerBoundaryForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexBilayerBoundaryForce)