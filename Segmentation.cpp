#include "Segmentation.h"\

/*!
 *  \brief Null constructor module.
 */
Segmentation::Segmentation()
{}

/*!
 *  \brief Constructor
 *  \param structure a reference to a ProteinStructure
 *  \param segments a reference to a vector<int>
 *  \param model_names a reference to a vector<string>
 */
Segmentation::Segmentation(ProteinStructure *structure, vector<int> &segments,
              vector<string> &model_names) : structure(structure), 
              segments(segments), model_names(model_names)
{}

