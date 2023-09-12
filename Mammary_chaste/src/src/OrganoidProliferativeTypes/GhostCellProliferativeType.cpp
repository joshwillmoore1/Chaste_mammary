#include "GhostCellProliferativeType.hpp"

GhostCellProliferativeType::GhostCellProliferativeType()
    : AbstractCellProliferativeType(200)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(GhostCellProliferativeType)
