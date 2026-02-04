#include "parameter_variant.h"
#include "../exception.h"

namespace ofec {
	//ParameterVariantStream& operator<< (ParameterVariantStream& stream, const ParameterVariant& var) {
	//	stream.input(var);
	//	return stream;
	//}
    std::string toString(ParameterType type) {
        switch (type) {
        case ParameterType::kBool:          return "kBool";
        case ParameterType::kInt:           return "kInt";
        case ParameterType::kSizeT:         return "kSizeT";
        case ParameterType::kChar:          return "kChar";
        case ParameterType::kReal:          return "kReal";
        case ParameterType::kVectorBool:    return "kVectorBool";
        case ParameterType::kVectorInt:     return "kVectorInt";
        case ParameterType::kVectorSizeT:   return "kVectorSizeT";
        case ParameterType::kString:        return "kString";
        case ParameterType::kVectorReal:    return "kVectorReal";
        default:                            return "Unknown";
        }
    }



    ParameterType fromString(const std::string& s) {
        if (s == "kBool")        return ParameterType::kBool;
        if (s == "kInt")         return ParameterType::kInt;
        if (s == "kSizeT")       return ParameterType::kSizeT;
        if (s == "kChar")        return ParameterType::kChar;
        if (s == "kReal")        return ParameterType::kReal;
        if (s == "kVectorBool")  return ParameterType::kVectorBool;
        if (s == "kVectorInt")   return ParameterType::kVectorInt;
        if (s == "kVectorSizeT") return ParameterType::kVectorSizeT;
        if (s == "kString")      return ParameterType::kString;
        if (s == "kVectorReal")  return ParameterType::kVectorReal;

        throw Exception("Unknown ParameterType: " + s);
    }

}