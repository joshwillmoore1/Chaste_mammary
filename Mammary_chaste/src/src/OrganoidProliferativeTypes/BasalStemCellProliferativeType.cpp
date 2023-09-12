#include "BasalStemCellProliferativeType.hpp"

BasalStemCellProliferativeType::BasalStemCellProliferativeType()
    : AbstractCellProliferativeType(10)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(BasalStemCellProliferativeType)
