#include "LumenERPositiveCellProliferativeType.hpp"

LumenERPositiveCellProliferativeType::LumenERPositiveCellProliferativeType()
    : AbstractCellProliferativeType(20)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(LumenERPositiveCellProliferativeType)
