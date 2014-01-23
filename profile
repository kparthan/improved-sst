Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls   s/call   s/call  name    
 14.76      2.00     2.00 137082418     0.00     0.00  convertToCartesian(double, double, double)
 11.25      3.53     1.53 26568565     0.00     0.00  convertToSpherical(lcb::Point<double>&)
  7.45      4.54     1.01 54433446     0.00     0.00  lcb::Vector<double>::l2Norm() const
  7.20      5.51     0.98 41151828     0.00     0.00  std::vector<double, std::allocator<double> >::vector(std::vector<double, std::allocator<double> > const&)
  6.05      6.33     0.82   324006     0.00     0.00  Superpose3DClass::diagonalize()
  6.01      7.15     0.82  2592056     0.00     0.00  void std::vector<Component, std::allocator<Component> >::_M_emplace_back_aux<Component const&>(Component const&)
  5.54      7.90     0.75 26568565     0.00     0.00  lcb::Vector<double>::normalize()
  5.35      8.62     0.73 25920480     0.00     0.00  Component::conflate(Component&)
  3.58      9.11     0.49   324006     0.00     0.00  Mixture::conflate(Component&)
  3.47      9.58     0.47   324006     0.00     0.00  Mixture::Mixture(int, std::vector<Component, std::allocator<Component> >&, std::vector<double, std::allocator<double> >&)
  3.47     10.05     0.47       16     0.03     0.03  std::_Sp_counted_base<(__gnu_cxx::_Lock_policy)2>::_M_release()
  2.66     10.41     0.36 26244566     0.00     0.00  VonMises3D::VonMises3D(std::array<double, 2ul>&, double)
  1.85     10.66     0.25 32428320     0.00     0.00  VonMises3D::density(double, double)
  1.70     10.89     0.23    15318     0.00     0.00  Segment::fitIdealModel(IdealModel&, Mixture&, int)
  1.62     11.11     0.22  4536960     0.00     0.00  lcb::Matrix<double>::operator*(lcb::Vector<double> const&) const
  1.18     11.27     0.16 26244566     0.00     0.00  VonMises3D::VonMises3D()
  1.07     11.41     0.15                             lcb::GenericData::~GenericData()
  1.03     11.55     0.14   405354     0.00     0.00  Mixture::probability(std::array<double, 2ul>&)
  1.03     11.69     0.14 11031522     0.00     0.00  std::vector<double, std::allocator<double> >::operator=(std::vector<double, std::allocator<double> > const&)
  0.89     11.81     0.12   704753     0.00     0.00  std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::operator=(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > > const&)
  0.81     11.92     0.11   648158     0.00     0.00  lcb::Matrix<double>::operator*(lcb::Matrix<double> const&)
  0.77     12.03     0.11   324079     0.00     0.00  convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&)
  0.59     12.11     0.08 26553254     0.00     0.00  VonMises3D::operator=(VonMises3D const&)
  0.59     12.19     0.08  1053396     0.00     0.00  void std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::_M_emplace_back_aux<std::vector<double, std::allocator<double> > const&>(std::vector<double, std::allocator<double> > const&)
  0.59     12.27     0.08   308688     0.00     0.00  Superpose3DClass::Superpose3DClass(suffStatClass&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&)
  0.55     12.34     0.08  3888948     0.00     0.00  lcb::Point<double> lcb::geometry::transform<double>(lcb::Point<double> const&, lcb::Matrix<double> const&)
  0.55     12.42     0.08  1944328     0.00     0.00  std::vector<double, std::allocator<double> >::_M_fill_insert(__gnu_cxx::__normal_iterator<double*, std::vector<double, std::allocator<double> > >, unsigned long, double const&)
  0.48     12.48     0.07  1620249     0.00     0.00  lcb::Matrix<double>::Matrix(lcb::Matrix<double> const&)
  0.44     12.54     0.06   972237     0.00     0.00  lcb::Matrix<double>::Matrix(int)
  0.44     12.60     0.06    15325     0.00     0.00  IdealModel::IdealModel(lcb::ProteinStructure*, std::string)
  0.44     12.66     0.06   648012     0.00     0.00  lcb::Matrix<double>::operator=(lcb::Matrix<double> const&)
  0.44     12.72     0.06   308688     0.00     0.00  Superpose3DClass::Superpose3DClass(std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&)
  0.41     12.78     0.06                             std::vector<lcb::Vector<double>, std::allocator<lcb::Vector<double> > >::~vector()
  0.37     12.83     0.05 32428320     0.00     0.00  Component::likelihood(std::array<double, 2ul>&)
  0.37     12.88     0.05   324006     0.00     0.00  Superpose3DClass::transformVectors(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&)
  0.30     12.92     0.04  2414275     0.00     0.00  void std::vector<double, std::allocator<double> >::_M_emplace_back_aux<double const&>(double const&)
  0.30     12.96     0.04        1     0.04     0.04  lcb::Chain::~Chain()
  0.30     13.00     0.04                             std::vector<lcb::Atom, std::allocator<lcb::Atom> >::~vector()
  0.26     13.03     0.04  1296316     0.00     0.00  std::vector<double, std::allocator<double> >::_M_default_append(unsigned long)
  0.22     13.06     0.03   648158     0.00     0.00  lcb::Vector<double>::angleBetween(lcb::Vector<double> const&, lcb::Vector<double> const&)
  0.22     13.09     0.03   354715     0.00     0.00  std::vector<lcb::Vector<double>, std::allocator<lcb::Vector<double> > >::_M_default_append(unsigned long)
  0.22     13.12     0.03   324006     0.00     0.00  Segment::getCurrentMeanAndDirection(std::pair<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >, lcb::Matrix<double> >&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&, int)
  0.22     13.15     0.03    26105     0.00     0.00  OptimalFit::operator=(OptimalFit const&)
  0.15     13.17     0.02  2268553     0.00     0.00  lcb::Matrix<double>::~Matrix()
  0.15     13.19     0.02   648158     0.00     0.00  lcb::Vector<double>::operator*(lcb::Vector<double> const&) const
  0.15     13.21     0.02   405354     0.00     0.00  std::vector<double, std::allocator<double> >::vector(unsigned long, double const&, std::allocator<double> const&)
  0.15     13.23     0.02   308688     0.00     0.00  Superpose3DClass::updateCenterOfMasses()
  0.15     13.25     0.02    17790     0.00     0.00  OptimalFit::OptimalFit(IdealModel&, double)
  0.15     13.27     0.02        1     0.02     0.02  extractName(std::string&)
  0.15     13.29     0.02                             Normal::negativeLogLikelihood(std::vector<double, std::allocator<double> >&)
  0.11     13.31     0.02   648158     0.00     0.00  lcb::Vector<double>::crossProduct(lcb::Vector<double> const&, lcb::Vector<double> const&)
  0.11     13.32     0.02   648158     0.00     0.00  lcb::Vector<double>::normalize_copy() const
  0.11     13.34     0.02                             assignSecondaryStructure(std::string, std::string, int)
  0.07     13.35     0.01   648158     0.00     0.00  lcb::Matrix<double> lcb::geometry::rotationMatrix<double>(lcb::Vector<double> const&, double)
  0.07     13.36     0.01   648012     0.00     0.00  lcb::Point<double> lcb::geometry::transform<double>(lcb::Point<double> const&, lcb::Matrix<double> const&)
  0.07     13.37     0.01   354715     0.00     0.00  void std::__uninitialized_default_n_1<false>::__uninit_default_n<lcb::Vector<double>*, unsigned long>(lcb::Vector<double>*, unsigned long)
  0.07     13.38     0.01   324006     0.00     0.00  Superpose3DClass::getSufficientStatistics()
  0.07     13.39     0.01   308688     0.00     0.00  Component::operator=(Component const&)
  0.07     13.40     0.01    26105     0.00     0.00  IdealModel::operator=(IdealModel const&)
  0.07     13.41     0.01    17794     0.00     0.00  Message::encodeUsingLogStarModel(double)
  0.07     13.42     0.01    15318     0.00     0.00  Superpose3DClass::computeQuaternionMatrix()
  0.07     13.43     0.01     2472     0.00     0.00  Segment::Segment(int, int, std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&, std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > >&)
  0.07     13.44     0.01       17     0.00     0.00  std::_Rb_tree<std::string, std::pair<std::string const, boost::program_options::variable_value>, std::_Select1st<std::pair<std::string const, boost::program_options::variable_value> >, std::less<std::string>, std::allocator<std::pair<std::string const, boost::program_options::variable_value> > >::find(std::string const&) const
  0.07     13.45     0.01        1     0.01    12.64  Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int)
  0.07     13.46     0.01        1     0.01     0.01  std::vector<IdealModel, std::allocator<IdealModel> >::~vector()
  0.07     13.47     0.01                             Superpose3DClass::getRMSD()
  0.07     13.48     0.01                             Component::computeSecondDerivative(double)
  0.07     13.49     0.01                             std::_Sp_counted_base<(__gnu_cxx::_Lock_policy)2>::~_Sp_counted_base()
  0.07     13.50     0.01                             std::vector<std::vector<OptimalFit, std::allocator<OptimalFit> >, std::allocator<std::vector<OptimalFit, std::allocator<OptimalFit> > > >::~vector()
  0.04     13.50     0.01   804749     0.00     0.00  point2vector(lcb::Point<double>&)
  0.04     13.51     0.01   648158     0.00     0.00  lcb::Matrix<double>::Matrix(int, int)
  0.04     13.51     0.01        8     0.00     0.00  void std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > >::_M_emplace_back_aux<std::array<double, 3ul> const&>(std::array<double, 3ul> const&)
  0.04     13.52     0.01        1     0.01     0.01  boost::program_options::basic_command_line_parser<char>::run()
  0.04     13.52     0.01                             generateRandomWeights(int, double)
  0.04     13.53     0.01                             lcb::Model::~Model()
  0.04     13.53     0.01                             boost::program_options::error_with_option_name::set_option_name(std::string const&)
  0.04     13.54     0.01                             Mixture::initialize2()
  0.04     13.54     0.01                             Component::generate(int)
  0.04     13.55     0.01                             std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >::~vector()
  0.04     13.55     0.01                             std::vector<double, std::allocator<double> >::resize(unsigned long)
  0.00     13.55     0.00   405354     0.00     0.00  Message::encodeUsingMixtureModel(std::array<double, 2ul>&, Mixture&)
  0.00     13.55     0.00   324086     0.00     0.00  Component::Component(std::array<double, 2ul>&, double, int)
  0.00     13.55     0.00    82139     0.00     0.00  Normal::density(double)
  0.00     13.55     0.00    81348     0.00     0.00  Message::encodeUsingNormalModel(double, Normal&)
  0.00     13.55     0.00    30636     0.00     0.00  lcb::Matrix<double>::Matrix()
  0.00     13.55     0.00    17792     0.00     0.00  Normal::Normal(double, double)
  0.00     13.55     0.00    17792     0.00     0.00  Message::Message()
  0.00     13.55     0.00    15319     0.00     0.00  Mixture::Mixture()
  0.00     13.55     0.00    15318     0.00     0.00  IdealModel::getResidues(int)
  0.00     13.55     0.00    15318     0.00     0.00  IdealModel::getStructure()
  0.00     13.55     0.00    15318     0.00     0.00  IdealModel::getName()
  0.00     13.55     0.00    15318     0.00     0.00  IdealModel::setLength(int)
  0.00     13.55     0.00    15318     0.00     0.00  OptimalFit::operator<(OptimalFit const&)
  0.00     13.55     0.00    15318     0.00     0.00  Superpose3DClass::Superpose3DClass(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&, std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&)
  0.00     13.55     0.00     5852     0.00     0.00  OptimalFit::OptimalFit(OptimalFit const&)
  0.00     13.55     0.00     4945     0.00     0.00  IdealModel::IdealModel()
  0.00     13.55     0.00     4945     0.00     0.00  OptimalFit::OptimalFit()
  0.00     13.55     0.00     2472     0.00     0.00  IdealModel::IdealModel(int, std::string)
  0.00     13.55     0.00     2472     0.00     0.00  Segment::fitNullModel(Mixture&)
  0.00     13.55     0.00     2472     0.00     0.00  OptimalFit::getMessageLength() const
  0.00     13.55     0.00      791     0.00     0.00  Message::encodeUsingSphereModel(double, Normal&)
  0.00     13.55     0.00      480     0.00     0.00  bool boost::char_separator<char, std::char_traits<char> >::operator()<__gnu_cxx::__normal_iterator<char const*, std::string>, std::string>(__gnu_cxx::__normal_iterator<char const*, std::string>&, __gnu_cxx::__normal_iterator<char const*, std::string>, std::string&)
  0.00     13.55     0.00      320     0.00     0.00  boost::char_separator<char, std::char_traits<char> >::~char_separator()
  0.00     13.55     0.00      320     0.00     0.00  boost::token_iterator<boost::char_separator<char, std::char_traits<char> >, __gnu_cxx::__normal_iterator<char const*, std::string>, std::string>::~token_iterator()
  0.00     13.55     0.00      240     0.00     0.00  boost::char_separator<char, std::char_traits<char> >::char_separator(boost::char_separator<char, std::char_traits<char> > const&)
  0.00     13.55     0.00      160     0.00     0.00  std::vector<double, std::allocator<double> >::push_back(double const&)
  0.00     13.55     0.00      120     0.00     0.00  std::_Rb_tree_node<std::pair<std::string const, std::string> >* std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_create_node<std::pair<std::string const, std::string> const&>(std::pair<std::string const, std::string> const&)
  0.00     13.55     0.00       80     0.00     0.00  std::vector<Component, std::allocator<Component> >::push_back(Component const&)
  0.00     13.55     0.00       76     0.00     0.00  int minimum<int>(int, int)
  0.00     13.55     0.00       73     0.00     0.00  Protein::computeTransformation(int, int)
  0.00     13.55     0.00       48     0.00     0.00  Segment::setInitialDistances(double, double)
  0.00     13.55     0.00       27     0.00     0.00  boost::detail::sp_counted_base::release()
  0.00     13.55     0.00       16     0.00     0.00  lcb::Model::~Model()
  0.00     13.55     0.00        9     0.00     0.00  std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, std::string> >*)
  0.00     13.55     0.00        8     0.00     0.06  parsePDBFile(std::string&)
  0.00     13.55     0.00        8     0.00     0.00  checkFile(std::string&)
  0.00     13.55     0.00        8     0.00     0.00  IdealModel::~IdealModel()
  0.00     13.55     0.00        8     0.00     0.00  boost::program_options::typed_value<std::string, char>* boost::program_options::value<std::string>(std::string*)
  0.00     13.55     0.00        8     0.00     0.00  lcb::GenericData::getIdentifier() const
  0.00     13.55     0.00        8     0.00     0.00  std::_Sp_counted_ptr_inplace<lcb::Model, std::allocator<lcb::Model>, (__gnu_cxx::_Lock_policy)2>::_M_dispose()
  0.00     13.55     0.00        8     0.00     0.00  std::_Sp_counted_ptr_inplace<lcb::Model, std::allocator<lcb::Model>, (__gnu_cxx::_Lock_policy)2>::_M_get_deleter(std::type_info const&)
  0.00     13.55     0.00        8     0.00     0.00  void std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >::_M_emplace_back_aux<lcb::Point<double> const&>(lcb::Point<double> const&)
  0.00     13.55     0.00        8     0.00     0.00  void std::vector<std::vector<OptimalFit, std::allocator<OptimalFit> >, std::allocator<std::vector<OptimalFit, std::allocator<OptimalFit> > > >::_M_emplace_back_aux<std::vector<OptimalFit, std::allocator<OptimalFit> > const&>(std::vector<OptimalFit, std::allocator<OptimalFit> > const&)
  0.00     13.55     0.00        8     0.00     0.00  void std::vector<std::string, std::allocator<std::string> >::_M_emplace_back_aux<std::string>(std::string&&)
  0.00     13.55     0.00        8     0.00     0.00  std::_Rb_tree<std::string, std::pair<std::string const, char>, std::_Select1st<std::pair<std::string const, char> >, std::less<std::string>, std::allocator<std::pair<std::string const, char> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, char> >*)
  0.00     13.55     0.00        8     0.00     0.00  std::_Rb_tree<std::string, std::pair<std::string const, int>, std::_Select1st<std::pair<std::string const, int> >, std::less<std::string>, std::allocator<std::pair<std::string const, int> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, int> >*)
  0.00     13.55     0.00        7     0.00     0.00  std::vector<IdealModel, std::allocator<IdealModel> >::push_back(IdealModel const&)
  0.00     13.55     0.00        5     0.00     0.00  std::vector<std::string, std::allocator<std::string> >::~vector()
  0.00     13.55     0.00        4     0.00     0.00  IdealModel::IdealModel(IdealModel const&)
  0.00     13.55     0.00        4     0.00     0.00  boost::program_options::typed_value<int, char>* boost::program_options::value<int>(int*)
  0.00     13.55     0.00        4     0.00     0.00  boost::function1<std::pair<std::string, std::string>, std::string const&>::clear()
  0.00     13.55     0.00        4     0.00     0.00  void std::vector<IdealModel, std::allocator<IdealModel> >::_M_emplace_back_aux<IdealModel const&>(IdealModel const&)
  0.00     13.55     0.00        4     0.00     0.00  std::vector<boost::program_options::basic_option<char>, std::allocator<boost::program_options::basic_option<char> > >::~vector()
  0.00     13.55     0.00        4     0.00     0.00  void std::vector<int, std::allocator<int> >::_M_emplace_back_aux<int const&>(int const&)
  0.00     13.55     0.00        2     0.00     0.00  writeToFile(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&, char const*)
  0.00     13.55     0.00        2     0.00     0.00  getPDBFilePath(std::string&)
  0.00     13.55     0.00        2     0.00     0.00  std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::~vector()
  0.00     13.55     0.00        2     0.00     0.00  void std::vector<int, std::allocator<int> >::_M_emplace_back_aux<int>(int&&)
  0.00     13.55     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN10IdealModelC2Ev
  0.00     13.55     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN10OptimalFitC2Ev
  0.00     13.55     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN10VonMises3DC2Ev
  0.00     13.55     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN16Superpose3DClass23computeRotationalCenterEv
  0.00     13.55     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN6NormalC2Ev
  0.00     13.55     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN7MessageC2Ev
  0.00     13.55     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN7MixtureC2Ev
  0.00     13.55     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN7ProteinC2Ev
  0.00     13.55     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN7SegmentC2EiiRSt6vectorIN3lcb5PointIdEESaIS3_EERS0_ISt5arrayIdLm3EESaIS8_EE
  0.00     13.55     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN9ComponentC2Ev
  0.00     13.55     0.00        1     0.00     0.00  _GLOBAL__sub_I_initialize_components_from_file
  0.00     13.55     0.00        1     0.00     0.00  _GLOBAL__sub_I_main
  0.00     13.55     0.00        1     0.00     0.41  loadIdealModels()
  0.00     13.55     0.00        1     0.00     0.00  std::ostream& lcb::operator<< <double>(std::ostream&, lcb::Point<double> const&)
  0.00     13.55     0.00        1     0.00     0.01  boost::program_options::basic_parsed_options<char> boost::program_options::parse_command_line<char>(int, char const* const*, boost::program_options::options_description const&, int, boost::function1<std::pair<std::string, std::string>, std::string const&>)
  0.00     13.55     0.00        1     0.00     0.00  char* boost::detail::lcast_put_unsigned<std::char_traits<char>, unsigned int, char>(unsigned int, char*)
  0.00     13.55     0.00        1     0.00     0.00  Mixture::load(std::string&)
  0.00     13.55     0.00        1     0.00     0.00  Mixture::~Mixture()
  0.00     13.55     0.00        1     0.00     0.00  Protein::checkChainBreak(std::string&, std::vector<lcb::Atom, std::allocator<lcb::Atom> >&)
  0.00     13.55     0.00        1     0.00     0.00  Protein::printCodeLengthMatrix(int)
  0.00     13.55     0.00        1     0.00    13.06  Protein::compressUsingIdealModels(Mixture&, int)
  0.00     13.55     0.00        1     0.00     0.00  Protein::translateProteinToOrigin(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&)
  0.00     13.55     0.00        1     0.00     0.00  Protein::computeOptimalSegmentation(int)
  0.00     13.55     0.00        1     0.00     0.00  Protein::computeSuccessiveDistances()
  0.00     13.55     0.00        1     0.00     0.00  Protein::initializeCodeLengthMatrices(int)
  0.00     13.55     0.00        1     0.00     0.01  Protein::computeSphericalTransformation()
  0.00     13.55     0.00        1     0.00     0.00  Protein::computeMessageLengthUsingNullModel(Mixture&)
  0.00     13.55     0.00        1     0.00     0.00  Protein::computeMessageLengthUsingSphereModel()
  0.00     13.55     0.00        1     0.00     0.04  Protein::Protein(lcb::ProteinStructure*, std::string&)
  0.00     13.55     0.00        1     0.00     0.00  Protein::~Protein()
  0.00     13.55     0.00        1     0.00     0.00  std::vector<OptimalFit, std::allocator<OptimalFit> >::~vector()
  0.00     13.55     0.00        1     0.00     0.00  std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >::vector(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > const&)
  0.00     13.55     0.00        1     0.00     0.00  void std::vector<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >, std::allocator<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > > >::_M_emplace_back_aux<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > const&>(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > const&)
  0.00     13.55     0.00        1     0.00     0.00  void std::vector<std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > >, std::allocator<std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > > > >::_M_emplace_back_aux<std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > > const&>(std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > > const&)
  0.00     13.55     0.00        1     0.00     0.00  void std::vector<std::string, std::allocator<std::string> >::_M_emplace_back_aux<std::string const&>(std::string const&)
  0.00     13.55     0.00        1     0.00     0.00  std::vector<int, std::allocator<int> >::operator=(std::vector<int, std::allocator<int> > const&)
  0.00     13.55     0.00        1     0.00     0.00  std::_Rb_tree<std::string, std::string, std::_Identity<std::string>, std::less<std::string>, std::allocator<std::string> >::_M_erase(std::_Rb_tree_node<std::string>*)
  0.00     13.55     0.00        1     0.00     0.00  std::_Rb_tree<std::string, std::pair<std::string const, boost::program_options::variable_value>, std::_Select1st<std::pair<std::string const, boost::program_options::variable_value> >, std::less<std::string>, std::allocator<std::pair<std::string const, boost::program_options::variable_value> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, boost::program_options::variable_value> >*)
  0.00     13.55     0.00        1     0.00     0.00  std::basic_ostream<char, std::char_traits<char> >& std::operator<< <std::char_traits<char> >(std::basic_ostream<char, std::char_traits<char> >&, char const*) [clone .constprop.477]

 %         the percentage of the total running time of the
time       program used by this function.

cumulative a running sum of the number of seconds accounted
 seconds   for by this function and those listed above it.

 self      the number of seconds accounted for by this
seconds    function alone.  This is the major sort for this
           listing.

calls      the number of times this function was invoked, if
           this function is profiled, else blank.
 
 self      the average number of milliseconds spent in this
ms/call    function per call, if this function is profiled,
	   else blank.

 total     the average number of milliseconds spent in this
ms/call    function and its descendents per call, if this 
	   function is profiled, else blank.

name       the name of the function.  This is the minor sort
           for this listing. The index shows the location of
	   the function in the gprof listing. If the index is
	   in parenthesis it shows where it would appear in
	   the gprof listing if it were to be printed.

Copyright (C) 2012 Free Software Foundation, Inc.

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.

		     Call graph (explanation follows)


granularity: each sample hit covers 2 byte(s) for 0.07% of 13.55 seconds

index % time    self  children    called     name
                                                 <spontaneous>
[1]     97.4    0.02   13.19                 assignSecondaryStructure(std::string, std::string, int) [1]
                0.00   13.06       1/1           Protein::compressUsingIdealModels(Mixture&, int) [2]
                0.00    0.06       1/8           parsePDBFile(std::string&) [23]
                0.00    0.04       1/1           Protein::Protein(lcb::ProteinStructure*, std::string&) [54]
                0.02    0.00       1/1           extractName(std::string&) [63]
                0.00    0.01       1/1           Protein::computeSphericalTransformation() [81]
                0.00    0.00       1/1           Protein::computeMessageLengthUsingNullModel(Mixture&) [94]
                0.00    0.00       1/1           Mixture::load(std::string&) [95]
                0.00    0.00       1/1           Protein::computeMessageLengthUsingSphereModel() [96]
                0.00    0.00       1/1           Protein::computeSuccessiveDistances() [98]
                0.00    0.00       1/15319       Mixture::Mixture() [111]
                0.00    0.00       1/1           Mixture::~Mixture() [171]
                0.00    0.00       1/1           Protein::~Protein() [176]
-----------------------------------------------
                0.00   13.06       1/1           assignSecondaryStructure(std::string, std::string, int) [1]
[2]     96.4    0.00   13.06       1         Protein::compressUsingIdealModels(Mixture&, int) [2]
                0.01   12.63       1/1           Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int) [3]
                0.00    0.41       1/1           loadIdealModels() [24]
                0.01    0.00       1/1           std::vector<IdealModel, std::allocator<IdealModel> >::~vector() [74]
                0.00    0.00       1/1           Protein::computeOptimalSegmentation(int) [175]
-----------------------------------------------
                0.01   12.63       1/1           Protein::compressUsingIdealModels(Mixture&, int) [2]
[3]     93.3    0.01   12.63       1         Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int) [3]
                0.23   12.22   15318/15318       Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
                0.00    0.12    2472/2472        Segment::fitNullModel(Mixture&) [33]
                0.03    0.02   26105/26105       OptimalFit::operator=(OptimalFit const&) [49]
                0.01    0.00    2472/2472        Segment::Segment(int, int, std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&, std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > >&) [68]
                0.00    0.00       1/1           Protein::initializeCodeLengthMatrices(int) [97]
                0.00    0.00   15318/15318       OptimalFit::operator<(OptimalFit const&) [115]
                0.00    0.00    4944/4945        OptimalFit::OptimalFit() [118]
                0.00    0.00    2472/2472        OptimalFit::getMessageLength() const [120]
                0.00    0.00      76/76          int minimum<int>(int, int) [129]
                0.00    0.00      48/48          Segment::setInitialDistances(double, double) [130]
                0.00    0.00       1/1           Protein::printCodeLengthMatrix(int) [173]
-----------------------------------------------
                0.23   12.22   15318/15318       Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int) [3]
[4]     91.9    0.23   12.22   15318         Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
                0.49    7.91  324006/324006      Mixture::conflate(Component&) [5]
                0.10    1.08  324006/324079      convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
                0.08    0.91  308688/308688      Superpose3DClass::Superpose3DClass(suffStatClass&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&) [12]
                0.00    0.81  354024/405354      Message::encodeUsingMixtureModel(std::array<double, 2ul>&, Mixture&) [14]
                0.03    0.18  324006/324006      Segment::getCurrentMeanAndDirection(std::pair<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >, lcb::Matrix<double> >&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&, int) [27]
                0.11    0.07  648012/704753      std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::operator=(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > > const&) [28]
                0.06    0.08  648012/648012      lcb::Matrix<double>::operator=(lcb::Matrix<double> const&) [32]
                0.06    0.05   15318/15325       IdealModel::IdealModel(lcb::ProteinStructure*, std::string) [34]
                0.00    0.06   15318/15318       Superpose3DClass::Superpose3DClass(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&, std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&) [45]
                0.05    0.00  324006/324006      Superpose3DClass::transformVectors(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&) [48]
                0.02    0.00   15318/17790       OptimalFit::OptimalFit(IdealModel&, double) [62]
                0.02    0.00  225617/1053396     void std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::_M_emplace_back_aux<std::vector<double, std::allocator<double> > const&>(std::vector<double, std::allocator<double> > const&) [40]
                0.00    0.02  324006/324086      Component::Component(std::array<double, 2ul>&, double, int) [66]
                0.01    0.00  308688/308688      Component::operator=(Component const&) [69]
                0.01    0.00  324006/324006      Superpose3DClass::getSufficientStatistics() [71]
                0.01    0.00   15318/17794       Message::encodeUsingLogStarModel(double) [79]
                0.00    0.01   15318/15318       IdealModel::getResidues(int) [80]
                0.00    0.00  324006/137082418     convertToCartesian(double, double, double) [8]
                0.00    0.00   30636/30636       lcb::Matrix<double>::Matrix() [108]
                0.00    0.00   30018/81348       Message::encodeUsingNormalModel(double, Normal&) [107]
                0.00    0.00   15318/17792       Normal::Normal(double, double) [109]
                0.00    0.00   15318/17792       Message::Message() [110]
                0.00    0.00   15318/15319       Mixture::Mixture() [111]
                0.00    0.00   15318/15318       IdealModel::getStructure() [112]
                0.00    0.00   15318/15318       IdealModel::getName() [113]
                0.00    0.00   15318/15318       IdealModel::setLength(int) [114]
                0.00    0.00     618/791         Message::encodeUsingSphereModel(double, Normal&) [121]
-----------------------------------------------
                0.49    7.91  324006/324006      Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[5]     62.0    0.49    7.91  324006         Mixture::conflate(Component&) [5]
                0.73    5.90 25920480/25920480     Component::conflate(Component&) [6]
                0.82    0.00 2592048/2592056     void std::vector<Component, std::allocator<Component> >::_M_emplace_back_aux<Component const&>(Component const&) [17]
                0.47    0.00  324006/324006      Mixture::Mixture(int, std::vector<Component, std::allocator<Component> >&, std::vector<double, std::allocator<double> >&) [21]
-----------------------------------------------
                0.73    5.90 25920480/25920480     Mixture::conflate(Component&) [5]
[6]     48.9    0.73    5.90 25920480         Component::conflate(Component&) [6]
                1.49    2.31 25920480/26568565     convertToSpherical(lcb::Point<double>&) [7]
                0.76    0.00 51840960/137082418     convertToCartesian(double, double, double) [8]
                0.36    0.38 25920480/26244566     VonMises3D::VonMises3D(std::array<double, 2ul>&, double) [18]
                0.16    0.38 25920480/26244566     VonMises3D::VonMises3D() [20]
                0.08    0.00 25920480/26553254     VonMises3D::operator=(VonMises3D const&) [39]
-----------------------------------------------
                0.00    0.00      73/26568565     Protein::computeSphericalTransformation() [81]
                0.04    0.06  648012/26568565     Segment::getCurrentMeanAndDirection(std::pair<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >, lcb::Matrix<double> >&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&, int) [27]
                1.49    2.31 25920480/26568565     Component::conflate(Component&) [6]
[7]     28.7    1.53    2.37 26568565         convertToSpherical(lcb::Point<double>&) [7]
                0.75    0.49 26568565/26568565     lcb::Vector<double>::normalize() [9]
                0.63    0.00 26568565/41151828     std::vector<double, std::allocator<double> >::vector(std::vector<double, std::allocator<double> > const&) [13]
                0.49    0.00 26568565/54433446     lcb::Vector<double>::l2Norm() const [11]
-----------------------------------------------
                0.00    0.00  324006/137082418     Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
                0.38    0.00 26244566/137082418     VonMises3D::VonMises3D() [20]
                0.38    0.00 26244566/137082418     VonMises3D::VonMises3D(std::array<double, 2ul>&, double) [18]
                0.47    0.00 32428320/137082418     VonMises3D::density(double, double) [19]
                0.76    0.00 51840960/137082418     Component::conflate(Component&) [6]
[8]     14.8    2.00    0.00 137082418         convertToCartesian(double, double, double) [8]
-----------------------------------------------
                0.75    0.49 26568565/26568565     convertToSpherical(lcb::Point<double>&) [7]
[9]      9.2    0.75    0.49 26568565         lcb::Vector<double>::normalize() [9]
                0.49    0.00 26568565/54433446     lcb::Vector<double>::l2Norm() const [11]
-----------------------------------------------
                0.00    0.00      73/324079      Protein::computeTransformation(int, int) [93]
                0.10    1.08  324006/324079      Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[10]     8.7    0.11    1.08  324079         convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
                0.08    0.33 3888948/3888948     lcb::Point<double> lcb::geometry::transform<double>(lcb::Point<double> const&, lcb::Matrix<double> const&) [25]
                0.11    0.05  648158/648158      lcb::Matrix<double>::operator*(lcb::Matrix<double> const&) [29]
                0.15    0.00 6481580/41151828     std::vector<double, std::allocator<double> >::vector(std::vector<double, std::allocator<double> > const&) [13]
                0.01    0.10  648158/648158      lcb::Matrix<double> lcb::geometry::rotationMatrix<double>(lcb::Vector<double> const&, double) [35]
                0.03    0.03  648158/648158      lcb::Vector<double>::angleBetween(lcb::Vector<double> const&, lcb::Vector<double> const&) [42]
                0.03    0.01  324079/354715      std::vector<lcb::Vector<double>, std::allocator<lcb::Vector<double> > >::_M_default_append(unsigned long) [52]
                0.04    0.00 1296316/1296316     std::vector<double, std::allocator<double> >::_M_default_append(unsigned long) [56]
                0.02    0.02  648158/648158      lcb::Vector<double>::crossProduct(lcb::Vector<double> const&, lcb::Vector<double> const&) [57]
                0.03    0.00 2268553/11031522     std::vector<double, std::allocator<double> >::operator=(std::vector<double, std::allocator<double> > const&) [31]
                0.01    0.01  324079/1620249     lcb::Matrix<double>::Matrix(lcb::Matrix<double> const&) [37]
                0.02    0.00  324079/972237      lcb::Matrix<double>::Matrix(int) [43]
                0.01    0.00 1620395/2268553     lcb::Matrix<double>::~Matrix() [58]
-----------------------------------------------
                0.01    0.00  648158/54433446     lcb::Vector<double>::angleBetween(lcb::Vector<double> const&, lcb::Vector<double> const&) [42]
                0.01    0.00  648158/54433446     lcb::Vector<double>::normalize_copy() const [44]
                0.49    0.00 26568565/54433446     convertToSpherical(lcb::Point<double>&) [7]
                0.49    0.00 26568565/54433446     lcb::Vector<double>::normalize() [9]
[11]     7.5    1.01    0.00 54433446         lcb::Vector<double>::l2Norm() const [11]
-----------------------------------------------
                0.08    0.91  308688/308688      Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[12]     7.3    0.08    0.91  308688         Superpose3DClass::Superpose3DClass(suffStatClass&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&) [12]
                0.78    0.00  308688/324006      Superpose3DClass::diagonalize() [16]
                0.06    0.05  308688/308688      Superpose3DClass::Superpose3DClass(std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&) [36]
                0.02    0.00  308688/308688      Superpose3DClass::updateCenterOfMasses() [61]
-----------------------------------------------
                0.02    0.00  648158/41151828     lcb::Vector<double>::crossProduct(lcb::Vector<double> const&, lcb::Vector<double> const&) [57]
                0.03    0.00 1296316/41151828     lcb::Vector<double>::normalize_copy() const [44]
                0.04    0.00 1620249/41151828     lcb::Matrix<double>::Matrix(lcb::Matrix<double> const&) [37]
                0.11    0.00 4536960/41151828     lcb::Matrix<double>::operator*(lcb::Vector<double> const&) const [26]
                0.15    0.00 6481580/41151828     convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
                0.63    0.00 26568565/41151828     convertToSpherical(lcb::Point<double>&) [7]
[13]     7.2    0.98    0.00 41151828         std::vector<double, std::allocator<double> >::vector(std::vector<double, std::allocator<double> > const&) [13]
-----------------------------------------------
                0.00    0.00      73/405354      Protein::computeMessageLengthUsingNullModel(Mixture&) [94]
                0.00    0.12   51257/405354      Segment::fitNullModel(Mixture&) [33]
                0.00    0.81  354024/405354      Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[14]     6.9    0.00    0.93  405354         Message::encodeUsingMixtureModel(std::array<double, 2ul>&, Mixture&) [14]
                0.14    0.79  405354/405354      Mixture::probability(std::array<double, 2ul>&) [15]
-----------------------------------------------
                0.14    0.79  405354/405354      Message::encodeUsingMixtureModel(std::array<double, 2ul>&, Mixture&) [14]
[15]     6.9    0.14    0.79  405354         Mixture::probability(std::array<double, 2ul>&) [15]
                0.25    0.47 32428320/32428320     VonMises3D::density(double, double) [19]
                0.05    0.00 32428320/32428320     Component::likelihood(std::array<double, 2ul>&) [47]
                0.02    0.00  405354/405354      std::vector<double, std::allocator<double> >::vector(unsigned long, double const&, std::allocator<double> const&) [60]
-----------------------------------------------
                0.04    0.00   15318/324006      Superpose3DClass::Superpose3DClass(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&, std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&) [45]
                0.78    0.00  308688/324006      Superpose3DClass::Superpose3DClass(suffStatClass&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&) [12]
[16]     6.1    0.82    0.00  324006         Superpose3DClass::diagonalize() [16]
-----------------------------------------------
                0.00    0.00       8/2592056     Mixture::load(std::string&) [95]
                0.82    0.00 2592048/2592056     Mixture::conflate(Component&) [5]
[17]     6.0    0.82    0.00 2592056         void std::vector<Component, std::allocator<Component> >::_M_emplace_back_aux<Component const&>(Component const&) [17]
-----------------------------------------------
                0.00    0.00  324086/26244566     Component::Component(std::array<double, 2ul>&, double, int) [66]
                0.36    0.38 25920480/26244566     Component::conflate(Component&) [6]
[18]     5.5    0.36    0.38 26244566         VonMises3D::VonMises3D(std::array<double, 2ul>&, double) [18]
                0.38    0.00 26244566/137082418     convertToCartesian(double, double, double) [8]
-----------------------------------------------
                0.25    0.47 32428320/32428320     Mixture::probability(std::array<double, 2ul>&) [15]
[19]     5.3    0.25    0.47 32428320         VonMises3D::density(double, double) [19]
                0.47    0.00 32428320/137082418     convertToCartesian(double, double, double) [8]
-----------------------------------------------
                0.00    0.00  324086/26244566     Component::Component(std::array<double, 2ul>&, double, int) [66]
                0.16    0.38 25920480/26244566     Component::conflate(Component&) [6]
[20]     4.0    0.16    0.38 26244566         VonMises3D::VonMises3D() [20]
                0.38    0.00 26244566/137082418     convertToCartesian(double, double, double) [8]
-----------------------------------------------
                0.47    0.00  324006/324006      Mixture::conflate(Component&) [5]
[21]     3.5    0.47    0.00  324006         Mixture::Mixture(int, std::vector<Component, std::allocator<Component> >&, std::vector<double, std::allocator<double> >&) [21]
-----------------------------------------------
                0.47    0.00      16/16          parsePDBFile(std::string&) [23]
[22]     3.5    0.47    0.00      16         std::_Sp_counted_base<(__gnu_cxx::_Lock_policy)2>::_M_release() [22]
                0.00    0.00       8/16          lcb::Model::~Model() [132]
                0.00    0.00       8/8           std::_Sp_counted_ptr_inplace<lcb::Model, std::allocator<lcb::Model>, (__gnu_cxx::_Lock_policy)2>::_M_dispose() [138]
-----------------------------------------------
                0.00    0.06       1/8           assignSecondaryStructure(std::string, std::string, int) [1]
                0.00    0.41       7/8           loadIdealModels() [24]
[23]     3.5    0.00    0.47       8         parsePDBFile(std::string&) [23]
                0.47    0.00      16/16          std::_Sp_counted_base<(__gnu_cxx::_Lock_policy)2>::_M_release() [22]
                0.00    0.00       8/8           checkFile(std::string&) [134]
                0.00    0.00       8/8           lcb::GenericData::getIdentifier() const [137]
                0.00    0.00       8/8           std::_Sp_counted_ptr_inplace<lcb::Model, std::allocator<lcb::Model>, (__gnu_cxx::_Lock_policy)2>::_M_get_deleter(std::type_info const&) [139]
                0.00    0.00       8/16          lcb::Model::~Model() [132]
                0.00    0.00       8/8           std::_Rb_tree<std::string, std::pair<std::string const, int>, std::_Select1st<std::pair<std::string const, int> >, std::less<std::string>, std::allocator<std::pair<std::string const, int> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, int> >*) [144]
                0.00    0.00       8/9           std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, std::string> >*) [133]
                0.00    0.00       8/8           std::_Rb_tree<std::string, std::pair<std::string const, char>, std::_Select1st<std::pair<std::string const, char> >, std::less<std::string>, std::allocator<std::pair<std::string const, char> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, char> >*) [143]
-----------------------------------------------
                0.00    0.41       1/1           Protein::compressUsingIdealModels(Mixture&, int) [2]
[24]     3.0    0.00    0.41       1         loadIdealModels() [24]
                0.00    0.41       7/8           parsePDBFile(std::string&) [23]
                0.00    0.00       7/15325       IdealModel::IdealModel(lcb::ProteinStructure*, std::string) [34]
                0.00    0.00       7/7           std::vector<IdealModel, std::allocator<IdealModel> >::push_back(IdealModel const&) [145]
                0.00    0.00       7/8           IdealModel::~IdealModel() [135]
                0.00    0.00       4/4           void std::vector<IdealModel, std::allocator<IdealModel> >::_M_emplace_back_aux<IdealModel const&>(IdealModel const&) [150]
-----------------------------------------------
                0.08    0.33 3888948/3888948     convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
[25]     3.0    0.08    0.33 3888948         lcb::Point<double> lcb::geometry::transform<double>(lcb::Point<double> const&, lcb::Matrix<double> const&) [25]
                0.19    0.09 3888948/4536960     lcb::Matrix<double>::operator*(lcb::Vector<double> const&) const [26]
                0.05    0.00 1296316/1944328     std::vector<double, std::allocator<double> >::_M_fill_insert(__gnu_cxx::__normal_iterator<double*, std::vector<double, std::allocator<double> > >, unsigned long, double const&) [41]
-----------------------------------------------
                0.03    0.02  648012/4536960     lcb::Point<double> lcb::geometry::transform<double>(lcb::Point<double> const&, lcb::Matrix<double> const&) [38]
                0.19    0.09 3888948/4536960     lcb::Point<double> lcb::geometry::transform<double>(lcb::Point<double> const&, lcb::Matrix<double> const&) [25]
[26]     2.4    0.22    0.11 4536960         lcb::Matrix<double>::operator*(lcb::Vector<double> const&) const [26]
                0.11    0.00 4536960/41151828     std::vector<double, std::allocator<double> >::vector(std::vector<double, std::allocator<double> > const&) [13]
-----------------------------------------------
                0.03    0.18  324006/324006      Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[27]     1.5    0.03    0.18  324006         Segment::getCurrentMeanAndDirection(std::pair<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >, lcb::Matrix<double> >&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&, int) [27]
                0.04    0.06  648012/26568565     convertToSpherical(lcb::Point<double>&) [7]
                0.01    0.07  648012/648012      lcb::Point<double> lcb::geometry::transform<double>(lcb::Point<double> const&, lcb::Matrix<double> const&) [38]
-----------------------------------------------
                0.00    0.00   26105/704753      IdealModel::operator=(IdealModel const&) [65]
                0.01    0.00   30636/704753      Superpose3DClass::Superpose3DClass(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&, std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&) [45]
                0.11    0.07  648012/704753      Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[28]     1.5    0.12    0.08  704753         std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::operator=(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > > const&) [28]
                0.08    0.00 6170921/11031522     std::vector<double, std::allocator<double> >::operator=(std::vector<double, std::allocator<double> > const&) [31]
-----------------------------------------------
                0.11    0.05  648158/648158      convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
[29]     1.2    0.11    0.05  648158         lcb::Matrix<double>::operator*(lcb::Matrix<double> const&) [29]
                0.03    0.02  648158/1620249     lcb::Matrix<double>::Matrix(lcb::Matrix<double> const&) [37]
                0.01    0.00  648158/2268553     lcb::Matrix<double>::~Matrix() [58]
                0.01    0.00  648158/648158      lcb::Matrix<double>::Matrix(int, int) [82]
-----------------------------------------------
                                                 <spontaneous>
[30]     1.1    0.15    0.00                 lcb::GenericData::~GenericData() [30]
-----------------------------------------------
                0.03    0.00 2268553/11031522     convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
                0.03    0.00 2592048/11031522     lcb::Matrix<double>::operator=(lcb::Matrix<double> const&) [32]
                0.08    0.00 6170921/11031522     std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::operator=(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > > const&) [28]
[31]     1.0    0.14    0.00 11031522         std::vector<double, std::allocator<double> >::operator=(std::vector<double, std::allocator<double> > const&) [31]
-----------------------------------------------
                0.06    0.08  648012/648012      Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[32]     1.0    0.06    0.08  648012         lcb::Matrix<double>::operator=(lcb::Matrix<double> const&) [32]
                0.03    0.02  648012/1620249     lcb::Matrix<double>::Matrix(lcb::Matrix<double> const&) [37]
                0.03    0.00 2592048/11031522     std::vector<double, std::allocator<double> >::operator=(std::vector<double, std::allocator<double> > const&) [31]
                0.00    0.00   30636/354715      std::vector<lcb::Vector<double>, std::allocator<lcb::Vector<double> > >::_M_default_append(unsigned long) [52]
-----------------------------------------------
                0.00    0.12    2472/2472        Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int) [3]
[33]     0.9    0.00    0.12    2472         Segment::fitNullModel(Mixture&) [33]
                0.00    0.12   51257/405354      Message::encodeUsingMixtureModel(std::array<double, 2ul>&, Mixture&) [14]
                0.00    0.00    2472/17790       OptimalFit::OptimalFit(IdealModel&, double) [62]
                0.00    0.00    2472/17794       Message::encodeUsingLogStarModel(double) [79]
                0.00    0.00   51257/81348       Message::encodeUsingNormalModel(double, Normal&) [107]
                0.00    0.00    2472/17792       Normal::Normal(double, double) [109]
                0.00    0.00    2472/17792       Message::Message() [110]
                0.00    0.00    2472/2472        IdealModel::IdealModel(int, std::string) [119]
                0.00    0.00      96/791         Message::encodeUsingSphereModel(double, Normal&) [121]
-----------------------------------------------
                0.00    0.00       7/15325       loadIdealModels() [24]
                0.06    0.05   15318/15325       Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[34]     0.8    0.06    0.05   15325         IdealModel::IdealModel(lcb::ProteinStructure*, std::string) [34]
                0.00    0.04  750925/804749      point2vector(lcb::Point<double>&) [50]
                0.01    0.00  107275/1053396     void std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::_M_emplace_back_aux<std::vector<double, std::allocator<double> > const&>(std::vector<double, std::allocator<double> > const&) [40]
-----------------------------------------------
                0.01    0.10  648158/648158      convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
[35]     0.8    0.01    0.10  648158         lcb::Matrix<double> lcb::geometry::rotationMatrix<double>(lcb::Vector<double> const&, double) [35]
                0.02    0.04  648158/648158      lcb::Vector<double>::normalize_copy() const [44]
                0.04    0.00  648158/972237      lcb::Matrix<double>::Matrix(int) [43]
-----------------------------------------------
                0.06    0.05  308688/308688      Superpose3DClass::Superpose3DClass(suffStatClass&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&) [12]
[36]     0.8    0.06    0.05  308688         Superpose3DClass::Superpose3DClass(std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&) [36]
                0.05    0.00  617376/1053396     void std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::_M_emplace_back_aux<std::vector<double, std::allocator<double> > const&>(std::vector<double, std::allocator<double> > const&) [40]
-----------------------------------------------
                0.01    0.01  324079/1620249     convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
                0.03    0.02  648012/1620249     lcb::Matrix<double>::operator=(lcb::Matrix<double> const&) [32]
                0.03    0.02  648158/1620249     lcb::Matrix<double>::operator*(lcb::Matrix<double> const&) [29]
[37]     0.8    0.07    0.04 1620249         lcb::Matrix<double>::Matrix(lcb::Matrix<double> const&) [37]
                0.04    0.00 1620249/41151828     std::vector<double, std::allocator<double> >::vector(std::vector<double, std::allocator<double> > const&) [13]
-----------------------------------------------
                0.01    0.07  648012/648012      Segment::getCurrentMeanAndDirection(std::pair<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >, lcb::Matrix<double> >&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&, int) [27]
[38]     0.6    0.01    0.07  648012         lcb::Point<double> lcb::geometry::transform<double>(lcb::Point<double> const&, lcb::Matrix<double> const&) [38]
                0.03    0.02  648012/4536960     lcb::Matrix<double>::operator*(lcb::Vector<double> const&) const [26]
                0.02    0.00  648012/1944328     std::vector<double, std::allocator<double> >::_M_fill_insert(__gnu_cxx::__normal_iterator<double*, std::vector<double, std::allocator<double> > >, unsigned long, double const&) [41]
-----------------------------------------------
                0.00    0.00  308688/26553254     Component::operator=(Component const&) [69]
                0.00    0.00  324086/26553254     Component::Component(std::array<double, 2ul>&, double, int) [66]
                0.08    0.00 25920480/26553254     Component::conflate(Component&) [6]
[39]     0.6    0.08    0.00 26553254         VonMises3D::operator=(VonMises3D const&) [39]
-----------------------------------------------
                0.00    0.00       1/1053396     Protein::computeSuccessiveDistances() [98]
                0.00    0.00       8/1053396     Protein::initializeCodeLengthMatrices(int) [97]
                0.00    0.00   13623/1053396     Segment::Segment(int, int, std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&, std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > >&) [68]
                0.01    0.00   89496/1053396     IdealModel::getResidues(int) [80]
                0.01    0.00  107275/1053396     IdealModel::IdealModel(lcb::ProteinStructure*, std::string) [34]
                0.02    0.00  225617/1053396     Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
                0.05    0.00  617376/1053396     Superpose3DClass::Superpose3DClass(std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&) [36]
[40]     0.6    0.08    0.00 1053396         void std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::_M_emplace_back_aux<std::vector<double, std::allocator<double> > const&>(std::vector<double, std::allocator<double> > const&) [40]
-----------------------------------------------
                0.02    0.00  648012/1944328     lcb::Point<double> lcb::geometry::transform<double>(lcb::Point<double> const&, lcb::Matrix<double> const&) [38]
                0.05    0.00 1296316/1944328     lcb::Point<double> lcb::geometry::transform<double>(lcb::Point<double> const&, lcb::Matrix<double> const&) [25]
[41]     0.6    0.08    0.00 1944328         std::vector<double, std::allocator<double> >::_M_fill_insert(__gnu_cxx::__normal_iterator<double*, std::vector<double, std::allocator<double> > >, unsigned long, double const&) [41]
-----------------------------------------------
                0.03    0.03  648158/648158      convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
[42]     0.5    0.03    0.03  648158         lcb::Vector<double>::angleBetween(lcb::Vector<double> const&, lcb::Vector<double> const&) [42]
                0.02    0.00  648158/648158      lcb::Vector<double>::operator*(lcb::Vector<double> const&) const [59]
                0.01    0.00  648158/54433446     lcb::Vector<double>::l2Norm() const [11]
-----------------------------------------------
                0.02    0.00  324079/972237      convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
                0.04    0.00  648158/972237      lcb::Matrix<double> lcb::geometry::rotationMatrix<double>(lcb::Vector<double> const&, double) [35]
[43]     0.4    0.06    0.00  972237         lcb::Matrix<double>::Matrix(int) [43]
-----------------------------------------------
                0.02    0.04  648158/648158      lcb::Matrix<double> lcb::geometry::rotationMatrix<double>(lcb::Vector<double> const&, double) [35]
[44]     0.4    0.02    0.04  648158         lcb::Vector<double>::normalize_copy() const [44]
                0.03    0.00 1296316/41151828     std::vector<double, std::allocator<double> >::vector(std::vector<double, std::allocator<double> > const&) [13]
                0.01    0.00  648158/54433446     lcb::Vector<double>::l2Norm() const [11]
-----------------------------------------------
                0.00    0.06   15318/15318       Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[45]     0.4    0.00    0.06   15318         Superpose3DClass::Superpose3DClass(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&, std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&) [45]
                0.04    0.00   15318/324006      Superpose3DClass::diagonalize() [16]
                0.01    0.00   15318/15318       Superpose3DClass::computeQuaternionMatrix() [72]
                0.01    0.00   30636/704753      std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::operator=(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > > const&) [28]
-----------------------------------------------
                                                 <spontaneous>
[46]     0.4    0.06    0.00                 std::vector<lcb::Vector<double>, std::allocator<lcb::Vector<double> > >::~vector() [46]
-----------------------------------------------
                0.05    0.00 32428320/32428320     Mixture::probability(std::array<double, 2ul>&) [15]
[47]     0.4    0.05    0.00 32428320         Component::likelihood(std::array<double, 2ul>&) [47]
-----------------------------------------------
                0.05    0.00  324006/324006      Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[48]     0.4    0.05    0.00  324006         Superpose3DClass::transformVectors(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&) [48]
-----------------------------------------------
                0.03    0.02   26105/26105       Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int) [3]
[49]     0.3    0.03    0.02   26105         OptimalFit::operator=(OptimalFit const&) [49]
                0.01    0.01   26105/26105       IdealModel::operator=(IdealModel const&) [65]
-----------------------------------------------
                0.00    0.00   53824/804749      Segment::Segment(int, int, std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&, std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > >&) [68]
                0.00    0.04  750925/804749      IdealModel::IdealModel(lcb::ProteinStructure*, std::string) [34]
[50]     0.3    0.01    0.04  804749         point2vector(lcb::Point<double>&) [50]
                0.04    0.00 2414247/2414275     void std::vector<double, std::allocator<double> >::_M_emplace_back_aux<double const&>(double const&) [51]
-----------------------------------------------
                0.00    0.00       8/2414275     Protein::computeSuccessiveDistances() [98]
                0.00    0.00      20/2414275     Mixture::load(std::string&) [95]
                0.04    0.00 2414247/2414275     point2vector(lcb::Point<double>&) [50]
[51]     0.3    0.04    0.00 2414275         void std::vector<double, std::allocator<double> >::_M_emplace_back_aux<double const&>(double const&) [51]
-----------------------------------------------
                0.00    0.00   30636/354715      lcb::Matrix<double>::operator=(lcb::Matrix<double> const&) [32]
                0.03    0.01  324079/354715      convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
[52]     0.3    0.03    0.01  354715         std::vector<lcb::Vector<double>, std::allocator<lcb::Vector<double> > >::_M_default_append(unsigned long) [52]
                0.01    0.00  354715/354715      void std::__uninitialized_default_n_1<false>::__uninit_default_n<lcb::Vector<double>*, unsigned long>(lcb::Vector<double>*, unsigned long) [70]
-----------------------------------------------
                0.04    0.00       1/1           Protein::Protein(lcb::ProteinStructure*, std::string&) [54]
[53]     0.3    0.04    0.00       1         lcb::Chain::~Chain() [53]
-----------------------------------------------
                0.00    0.04       1/1           assignSecondaryStructure(std::string, std::string, int) [1]
[54]     0.3    0.00    0.04       1         Protein::Protein(lcb::ProteinStructure*, std::string&) [54]
                0.04    0.00       1/1           lcb::Chain::~Chain() [53]
                0.00    0.00       8/8           void std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >::_M_emplace_back_aux<lcb::Point<double> const&>(lcb::Point<double> const&) [140]
                0.00    0.00       1/1           Protein::checkChainBreak(std::string&, std::vector<lcb::Atom, std::allocator<lcb::Atom> >&) [172]
                0.00    0.00       1/2           writeToFile(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&, char const*) [153]
                0.00    0.00       1/1           Protein::translateProteinToOrigin(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [174]
                0.00    0.00       1/1           void std::vector<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >, std::allocator<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > > >::_M_emplace_back_aux<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > const&>(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > const&) [179]
                0.00    0.00       1/1           void std::vector<std::string, std::allocator<std::string> >::_M_emplace_back_aux<std::string const&>(std::string const&) [181]
                0.00    0.00       1/5           std::vector<std::string, std::allocator<std::string> >::~vector() [146]
-----------------------------------------------
                                                 <spontaneous>
[55]     0.3    0.04    0.00                 std::vector<lcb::Atom, std::allocator<lcb::Atom> >::~vector() [55]
-----------------------------------------------
                0.04    0.00 1296316/1296316     convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
[56]     0.3    0.04    0.00 1296316         std::vector<double, std::allocator<double> >::_M_default_append(unsigned long) [56]
-----------------------------------------------
                0.02    0.02  648158/648158      convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
[57]     0.2    0.02    0.02  648158         lcb::Vector<double>::crossProduct(lcb::Vector<double> const&, lcb::Vector<double> const&) [57]
                0.02    0.00  648158/41151828     std::vector<double, std::allocator<double> >::vector(std::vector<double, std::allocator<double> > const&) [13]
-----------------------------------------------
                0.01    0.00  648158/2268553     lcb::Matrix<double>::operator*(lcb::Matrix<double> const&) [29]
                0.01    0.00 1620395/2268553     convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
[58]     0.1    0.02    0.00 2268553         lcb::Matrix<double>::~Matrix() [58]
-----------------------------------------------
                0.02    0.00  648158/648158      lcb::Vector<double>::angleBetween(lcb::Vector<double> const&, lcb::Vector<double> const&) [42]
[59]     0.1    0.02    0.00  648158         lcb::Vector<double>::operator*(lcb::Vector<double> const&) const [59]
-----------------------------------------------
                0.02    0.00  405354/405354      Mixture::probability(std::array<double, 2ul>&) [15]
[60]     0.1    0.02    0.00  405354         std::vector<double, std::allocator<double> >::vector(unsigned long, double const&, std::allocator<double> const&) [60]
-----------------------------------------------
                0.02    0.00  308688/308688      Superpose3DClass::Superpose3DClass(suffStatClass&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&) [12]
[61]     0.1    0.02    0.00  308688         Superpose3DClass::updateCenterOfMasses() [61]
-----------------------------------------------
                0.00    0.00    2472/17790       Segment::fitNullModel(Mixture&) [33]
                0.02    0.00   15318/17790       Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[62]     0.1    0.02    0.00   17790         OptimalFit::OptimalFit(IdealModel&, double) [62]
-----------------------------------------------
                0.02    0.00       1/1           assignSecondaryStructure(std::string, std::string, int) [1]
[63]     0.1    0.02    0.00       1         extractName(std::string&) [63]
-----------------------------------------------
                                                 <spontaneous>
[64]     0.1    0.02    0.00                 Normal::negativeLogLikelihood(std::vector<double, std::allocator<double> >&) [64]
-----------------------------------------------
                0.01    0.01   26105/26105       OptimalFit::operator=(OptimalFit const&) [49]
[65]     0.1    0.01    0.01   26105         IdealModel::operator=(IdealModel const&) [65]
                0.00    0.00   26105/704753      std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::operator=(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > > const&) [28]
-----------------------------------------------
                0.00    0.00      80/324086      Mixture::load(std::string&) [95]
                0.00    0.02  324006/324086      Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[66]     0.1    0.00    0.02  324086         Component::Component(std::array<double, 2ul>&, double, int) [66]
                0.00    0.00  324086/26244566     VonMises3D::VonMises3D(std::array<double, 2ul>&, double) [18]
                0.00    0.00  324086/26244566     VonMises3D::VonMises3D() [20]
                0.00    0.00  324086/26553254     VonMises3D::operator=(VonMises3D const&) [39]
-----------------------------------------------
                                                 <spontaneous>
[67]     0.1    0.00    0.02                 parseCommandLineInput(int, char**) [67]
                0.01    0.00      17/17          std::_Rb_tree<std::string, std::pair<std::string const, boost::program_options::variable_value>, std::_Select1st<std::pair<std::string const, boost::program_options::variable_value> >, std::less<std::string>, std::allocator<std::pair<std::string const, boost::program_options::variable_value> > >::find(std::string const&) const [73]
                0.00    0.01       1/1           boost::program_options::basic_parsed_options<char> boost::program_options::parse_command_line<char>(int, char const* const*, boost::program_options::options_description const&, int, boost::function1<std::pair<std::string, std::string>, std::string const&>) [84]
                0.00    0.00      24/27          boost::detail::sp_counted_base::release() [131]
                0.00    0.00       8/8           boost::program_options::typed_value<std::string, char>* boost::program_options::value<std::string>(std::string*) [136]
                0.00    0.00       4/4           boost::program_options::typed_value<int, char>* boost::program_options::value<int>(int*) [148]
                0.00    0.00       2/2           getPDBFilePath(std::string&) [154]
                0.00    0.00       1/4           std::vector<boost::program_options::basic_option<char>, std::allocator<boost::program_options::basic_option<char> > >::~vector() [151]
                0.00    0.00       1/4           boost::function1<std::pair<std::string, std::string>, std::string const&>::clear() [149]
                0.00    0.00       1/9           std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, std::string> >*) [133]
                0.00    0.00       1/1           std::_Rb_tree<std::string, std::string, std::_Identity<std::string>, std::less<std::string>, std::allocator<std::string> >::_M_erase(std::_Rb_tree_node<std::string>*) [183]
                0.00    0.00       1/1           std::_Rb_tree<std::string, std::pair<std::string const, boost::program_options::variable_value>, std::_Select1st<std::pair<std::string const, boost::program_options::variable_value> >, std::less<std::string>, std::allocator<std::pair<std::string const, boost::program_options::variable_value> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, boost::program_options::variable_value> >*) [184]
                0.00    0.00       1/1           std::basic_ostream<char, std::char_traits<char> >& std::operator<< <std::char_traits<char> >(std::basic_ostream<char, std::char_traits<char> >&, char const*) [clone .constprop.477] [185]
-----------------------------------------------
                0.01    0.00    2472/2472        Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int) [3]
[68]     0.1    0.01    0.00    2472         Segment::Segment(int, int, std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&, std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > >&) [68]
                0.00    0.00   53824/804749      point2vector(lcb::Point<double>&) [50]
                0.00    0.00   13623/1053396     void std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::_M_emplace_back_aux<std::vector<double, std::allocator<double> > const&>(std::vector<double, std::allocator<double> > const&) [40]
-----------------------------------------------
                0.01    0.00  308688/308688      Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[69]     0.1    0.01    0.00  308688         Component::operator=(Component const&) [69]
                0.00    0.00  308688/26553254     VonMises3D::operator=(VonMises3D const&) [39]
-----------------------------------------------
                0.01    0.00  354715/354715      std::vector<lcb::Vector<double>, std::allocator<lcb::Vector<double> > >::_M_default_append(unsigned long) [52]
[70]     0.1    0.01    0.00  354715         void std::__uninitialized_default_n_1<false>::__uninit_default_n<lcb::Vector<double>*, unsigned long>(lcb::Vector<double>*, unsigned long) [70]
-----------------------------------------------
                0.01    0.00  324006/324006      Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[71]     0.1    0.01    0.00  324006         Superpose3DClass::getSufficientStatistics() [71]
-----------------------------------------------
                0.01    0.00   15318/15318       Superpose3DClass::Superpose3DClass(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&, std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&) [45]
[72]     0.1    0.01    0.00   15318         Superpose3DClass::computeQuaternionMatrix() [72]
-----------------------------------------------
                0.01    0.00      17/17          parseCommandLineInput(int, char**) [67]
[73]     0.1    0.01    0.00      17         std::_Rb_tree<std::string, std::pair<std::string const, boost::program_options::variable_value>, std::_Select1st<std::pair<std::string const, boost::program_options::variable_value> >, std::less<std::string>, std::allocator<std::pair<std::string const, boost::program_options::variable_value> > >::find(std::string const&) const [73]
-----------------------------------------------
                0.01    0.00       1/1           Protein::compressUsingIdealModels(Mixture&, int) [2]
[74]     0.1    0.01    0.00       1         std::vector<IdealModel, std::allocator<IdealModel> >::~vector() [74]
-----------------------------------------------
                                                 <spontaneous>
[75]     0.1    0.01    0.00                 Superpose3DClass::getRMSD() [75]
-----------------------------------------------
                                                 <spontaneous>
[76]     0.1    0.01    0.00                 Component::computeSecondDerivative(double) [76]
-----------------------------------------------
                                                 <spontaneous>
[77]     0.1    0.01    0.00                 std::_Sp_counted_base<(__gnu_cxx::_Lock_policy)2>::~_Sp_counted_base() [77]
-----------------------------------------------
                                                 <spontaneous>
[78]     0.1    0.01    0.00                 std::vector<std::vector<OptimalFit, std::allocator<OptimalFit> >, std::allocator<std::vector<OptimalFit, std::allocator<OptimalFit> > > >::~vector() [78]
-----------------------------------------------
                0.00    0.00       2/17794       Protein::computeMessageLengthUsingSphereModel() [96]
                0.00    0.00       2/17794       Protein::computeMessageLengthUsingNullModel(Mixture&) [94]
                0.00    0.00    2472/17794       Segment::fitNullModel(Mixture&) [33]
                0.01    0.00   15318/17794       Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[79]     0.1    0.01    0.00   17794         Message::encodeUsingLogStarModel(double) [79]
-----------------------------------------------
                0.00    0.01   15318/15318       Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[80]     0.1    0.00    0.01   15318         IdealModel::getResidues(int) [80]
                0.01    0.00   89496/1053396     void std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::_M_emplace_back_aux<std::vector<double, std::allocator<double> > const&>(std::vector<double, std::allocator<double> > const&) [40]
-----------------------------------------------
                0.00    0.01       1/1           assignSecondaryStructure(std::string, std::string, int) [1]
[81]     0.0    0.00    0.01       1         Protein::computeSphericalTransformation() [81]
                0.01    0.00       8/8           void std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > >::_M_emplace_back_aux<std::array<double, 3ul> const&>(std::array<double, 3ul> const&) [83]
                0.00    0.00      73/73          Protein::computeTransformation(int, int) [93]
                0.00    0.00      73/26568565     convertToSpherical(lcb::Point<double>&) [7]
                0.00    0.00       1/1           void std::vector<std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > >, std::allocator<std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > > > >::_M_emplace_back_aux<std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > > const&>(std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > > const&) [180]
-----------------------------------------------
                0.01    0.00  648158/648158      lcb::Matrix<double>::operator*(lcb::Matrix<double> const&) [29]
[82]     0.0    0.01    0.00  648158         lcb::Matrix<double>::Matrix(int, int) [82]
-----------------------------------------------
                0.01    0.00       8/8           Protein::computeSphericalTransformation() [81]
[83]     0.0    0.01    0.00       8         void std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > >::_M_emplace_back_aux<std::array<double, 3ul> const&>(std::array<double, 3ul> const&) [83]
-----------------------------------------------
                0.00    0.01       1/1           parseCommandLineInput(int, char**) [67]
[84]     0.0    0.00    0.01       1         boost::program_options::basic_parsed_options<char> boost::program_options::parse_command_line<char>(int, char const* const*, boost::program_options::options_description const&, int, boost::function1<std::pair<std::string, std::string>, std::string const&>) [84]
                0.01    0.00       1/1           boost::program_options::basic_command_line_parser<char>::run() [85]
                0.00    0.00       8/8           void std::vector<std::string, std::allocator<std::string> >::_M_emplace_back_aux<std::string>(std::string&&) [142]
                0.00    0.00       3/5           std::vector<std::string, std::allocator<std::string> >::~vector() [146]
                0.00    0.00       3/4           boost::function1<std::pair<std::string, std::string>, std::string const&>::clear() [149]
-----------------------------------------------
                0.01    0.00       1/1           boost::program_options::basic_parsed_options<char> boost::program_options::parse_command_line<char>(int, char const* const*, boost::program_options::options_description const&, int, boost::function1<std::pair<std::string, std::string>, std::string const&>) [84]
[85]     0.0    0.01    0.00       1         boost::program_options::basic_command_line_parser<char>::run() [85]
                0.00    0.00       3/4           std::vector<boost::program_options::basic_option<char>, std::allocator<boost::program_options::basic_option<char> > >::~vector() [151]
-----------------------------------------------
                                                 <spontaneous>
[86]     0.0    0.01    0.00                 generateRandomWeights(int, double) [86]
-----------------------------------------------
                                                 <spontaneous>
[87]     0.0    0.01    0.00                 lcb::Model::~Model() [87]
-----------------------------------------------
                                                 <spontaneous>
[88]     0.0    0.01    0.00                 boost::program_options::error_with_option_name::set_option_name(std::string const&) [88]
-----------------------------------------------
                                                 <spontaneous>
[89]     0.0    0.01    0.00                 Mixture::initialize2() [89]
-----------------------------------------------
                                                 <spontaneous>
[90]     0.0    0.01    0.00                 Component::generate(int) [90]
-----------------------------------------------
                                                 <spontaneous>
[91]     0.0    0.01    0.00                 std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >::~vector() [91]
-----------------------------------------------
                                                 <spontaneous>
[92]     0.0    0.01    0.00                 std::vector<double, std::allocator<double> >::resize(unsigned long) [92]
-----------------------------------------------
                0.00    0.00      73/73          Protein::computeSphericalTransformation() [81]
[93]     0.0    0.00    0.00      73         Protein::computeTransformation(int, int) [93]
                0.00    0.00      73/324079      convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [10]
-----------------------------------------------
                0.00    0.00       1/1           assignSecondaryStructure(std::string, std::string, int) [1]
[94]     0.0    0.00    0.00       1         Protein::computeMessageLengthUsingNullModel(Mixture&) [94]
                0.00    0.00      73/405354      Message::encodeUsingMixtureModel(std::array<double, 2ul>&, Mixture&) [14]
                0.00    0.00       2/17794       Message::encodeUsingLogStarModel(double) [79]
                0.00    0.00      73/81348       Message::encodeUsingNormalModel(double, Normal&) [107]
                0.00    0.00       2/791         Message::encodeUsingSphereModel(double, Normal&) [121]
                0.00    0.00       1/17792       Message::Message() [110]
                0.00    0.00       1/17792       Normal::Normal(double, double) [109]
-----------------------------------------------
                0.00    0.00       1/1           assignSecondaryStructure(std::string, std::string, int) [1]
[95]     0.0    0.00    0.00       1         Mixture::load(std::string&) [95]
                0.00    0.00      80/324086      Component::Component(std::array<double, 2ul>&, double, int) [66]
                0.00    0.00       8/2592056     void std::vector<Component, std::allocator<Component> >::_M_emplace_back_aux<Component const&>(Component const&) [17]
                0.00    0.00      20/2414275     void std::vector<double, std::allocator<double> >::_M_emplace_back_aux<double const&>(double const&) [51]
                0.00    0.00     480/480         bool boost::char_separator<char, std::char_traits<char> >::operator()<__gnu_cxx::__normal_iterator<char const*, std::string>, std::string>(__gnu_cxx::__normal_iterator<char const*, std::string>&, __gnu_cxx::__normal_iterator<char const*, std::string>, std::string&) [122]
                0.00    0.00     320/320         boost::char_separator<char, std::char_traits<char> >::~char_separator() [123]
                0.00    0.00     320/320         boost::token_iterator<boost::char_separator<char, std::char_traits<char> >, __gnu_cxx::__normal_iterator<char const*, std::string>, std::string>::~token_iterator() [124]
                0.00    0.00     240/240         boost::char_separator<char, std::char_traits<char> >::char_separator(boost::char_separator<char, std::char_traits<char> > const&) [125]
                0.00    0.00     160/160         std::vector<double, std::allocator<double> >::push_back(double const&) [126]
                0.00    0.00      80/80          std::vector<Component, std::allocator<Component> >::push_back(Component const&) [128]
-----------------------------------------------
                0.00    0.00       1/1           assignSecondaryStructure(std::string, std::string, int) [1]
[96]     0.0    0.00    0.00       1         Protein::computeMessageLengthUsingSphereModel() [96]
                0.00    0.00       2/17794       Message::encodeUsingLogStarModel(double) [79]
                0.00    0.00      75/791         Message::encodeUsingSphereModel(double, Normal&) [121]
                0.00    0.00       1/17792       Normal::Normal(double, double) [109]
                0.00    0.00       1/17792       Message::Message() [110]
-----------------------------------------------
                0.00    0.00       1/1           Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int) [3]
[97]     0.0    0.00    0.00       1         Protein::initializeCodeLengthMatrices(int) [97]
                0.00    0.00       8/1053396     void std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::_M_emplace_back_aux<std::vector<double, std::allocator<double> > const&>(std::vector<double, std::allocator<double> > const&) [40]
                0.00    0.00    5244/5852        OptimalFit::OptimalFit(OptimalFit const&) [116]
                0.00    0.00       8/8           void std::vector<std::vector<OptimalFit, std::allocator<OptimalFit> >, std::allocator<std::vector<OptimalFit, std::allocator<OptimalFit> > > >::_M_emplace_back_aux<std::vector<OptimalFit, std::allocator<OptimalFit> > const&>(std::vector<OptimalFit, std::allocator<OptimalFit> > const&) [141]
                0.00    0.00       1/4945        OptimalFit::OptimalFit() [118]
                0.00    0.00       1/1           std::vector<OptimalFit, std::allocator<OptimalFit> >::~vector() [177]
                0.00    0.00       1/8           IdealModel::~IdealModel() [135]
-----------------------------------------------
                0.00    0.00       1/1           assignSecondaryStructure(std::string, std::string, int) [1]
[98]     0.0    0.00    0.00       1         Protein::computeSuccessiveDistances() [98]
                0.00    0.00       8/2414275     void std::vector<double, std::allocator<double> >::_M_emplace_back_aux<double const&>(double const&) [51]
                0.00    0.00       1/1053396     void std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::_M_emplace_back_aux<std::vector<double, std::allocator<double> > const&>(std::vector<double, std::allocator<double> > const&) [40]
-----------------------------------------------
                0.00    0.00     791/82139       Message::encodeUsingSphereModel(double, Normal&) [121]
                0.00    0.00   81348/82139       Message::encodeUsingNormalModel(double, Normal&) [107]
[106]    0.0    0.00    0.00   82139         Normal::density(double) [106]
-----------------------------------------------
                0.00    0.00      73/81348       Protein::computeMessageLengthUsingNullModel(Mixture&) [94]
                0.00    0.00   30018/81348       Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
                0.00    0.00   51257/81348       Segment::fitNullModel(Mixture&) [33]
[107]    0.0    0.00    0.00   81348         Message::encodeUsingNormalModel(double, Normal&) [107]
                0.00    0.00   81348/82139       Normal::density(double) [106]
-----------------------------------------------
                0.00    0.00   30636/30636       Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[108]    0.0    0.00    0.00   30636         lcb::Matrix<double>::Matrix() [108]
-----------------------------------------------
                0.00    0.00       1/17792       Protein::computeMessageLengthUsingSphereModel() [96]
                0.00    0.00       1/17792       Protein::computeMessageLengthUsingNullModel(Mixture&) [94]
                0.00    0.00    2472/17792       Segment::fitNullModel(Mixture&) [33]
                0.00    0.00   15318/17792       Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[109]    0.0    0.00    0.00   17792         Normal::Normal(double, double) [109]
-----------------------------------------------
                0.00    0.00       1/17792       Protein::computeMessageLengthUsingSphereModel() [96]
                0.00    0.00       1/17792       Protein::computeMessageLengthUsingNullModel(Mixture&) [94]
                0.00    0.00    2472/17792       Segment::fitNullModel(Mixture&) [33]
                0.00    0.00   15318/17792       Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[110]    0.0    0.00    0.00   17792         Message::Message() [110]
-----------------------------------------------
                0.00    0.00       1/15319       assignSecondaryStructure(std::string, std::string, int) [1]
                0.00    0.00   15318/15319       Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[111]    0.0    0.00    0.00   15319         Mixture::Mixture() [111]
-----------------------------------------------
                0.00    0.00   15318/15318       Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[112]    0.0    0.00    0.00   15318         IdealModel::getStructure() [112]
-----------------------------------------------
                0.00    0.00   15318/15318       Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[113]    0.0    0.00    0.00   15318         IdealModel::getName() [113]
-----------------------------------------------
                0.00    0.00   15318/15318       Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[114]    0.0    0.00    0.00   15318         IdealModel::setLength(int) [114]
-----------------------------------------------
                0.00    0.00   15318/15318       Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int) [3]
[115]    0.0    0.00    0.00   15318         OptimalFit::operator<(OptimalFit const&) [115]
-----------------------------------------------
                0.00    0.00     608/5852        void std::vector<std::vector<OptimalFit, std::allocator<OptimalFit> >, std::allocator<std::vector<OptimalFit, std::allocator<OptimalFit> > > >::_M_emplace_back_aux<std::vector<OptimalFit, std::allocator<OptimalFit> > const&>(std::vector<OptimalFit, std::allocator<OptimalFit> > const&) [141]
                0.00    0.00    5244/5852        Protein::initializeCodeLengthMatrices(int) [97]
[116]    0.0    0.00    0.00    5852         OptimalFit::OptimalFit(OptimalFit const&) [116]
-----------------------------------------------
                0.00    0.00    4945/4945        OptimalFit::OptimalFit() [118]
[117]    0.0    0.00    0.00    4945         IdealModel::IdealModel() [117]
-----------------------------------------------
                0.00    0.00       1/4945        Protein::initializeCodeLengthMatrices(int) [97]
                0.00    0.00    4944/4945        Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int) [3]
[118]    0.0    0.00    0.00    4945         OptimalFit::OptimalFit() [118]
                0.00    0.00    4945/4945        IdealModel::IdealModel() [117]
-----------------------------------------------
                0.00    0.00    2472/2472        Segment::fitNullModel(Mixture&) [33]
[119]    0.0    0.00    0.00    2472         IdealModel::IdealModel(int, std::string) [119]
-----------------------------------------------
                0.00    0.00    2472/2472        Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int) [3]
[120]    0.0    0.00    0.00    2472         OptimalFit::getMessageLength() const [120]
-----------------------------------------------
                0.00    0.00       2/791         Protein::computeMessageLengthUsingNullModel(Mixture&) [94]
                0.00    0.00      75/791         Protein::computeMessageLengthUsingSphereModel() [96]
                0.00    0.00      96/791         Segment::fitNullModel(Mixture&) [33]
                0.00    0.00     618/791         Segment::fitIdealModel(IdealModel&, Mixture&, int) [4]
[121]    0.0    0.00    0.00     791         Message::encodeUsingSphereModel(double, Normal&) [121]
                0.00    0.00     791/82139       Normal::density(double) [106]
-----------------------------------------------
                0.00    0.00     480/480         Mixture::load(std::string&) [95]
[122]    0.0    0.00    0.00     480         bool boost::char_separator<char, std::char_traits<char> >::operator()<__gnu_cxx::__normal_iterator<char const*, std::string>, std::string>(__gnu_cxx::__normal_iterator<char const*, std::string>&, __gnu_cxx::__normal_iterator<char const*, std::string>, std::string&) [122]
-----------------------------------------------
                0.00    0.00     320/320         Mixture::load(std::string&) [95]
[123]    0.0    0.00    0.00     320         boost::char_separator<char, std::char_traits<char> >::~char_separator() [123]
-----------------------------------------------
                0.00    0.00     320/320         Mixture::load(std::string&) [95]
[124]    0.0    0.00    0.00     320         boost::token_iterator<boost::char_separator<char, std::char_traits<char> >, __gnu_cxx::__normal_iterator<char const*, std::string>, std::string>::~token_iterator() [124]
-----------------------------------------------
                0.00    0.00     240/240         Mixture::load(std::string&) [95]
[125]    0.0    0.00    0.00     240         boost::char_separator<char, std::char_traits<char> >::char_separator(boost::char_separator<char, std::char_traits<char> > const&) [125]
-----------------------------------------------
                0.00    0.00     160/160         Mixture::load(std::string&) [95]
[126]    0.0    0.00    0.00     160         std::vector<double, std::allocator<double> >::push_back(double const&) [126]
-----------------------------------------------
                0.00    0.00     120/120         std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_copy(std::_Rb_tree_node<std::pair<std::string const, std::string> > const*, std::_Rb_tree_node<std::pair<std::string const, std::string> >*) [580]
[127]    0.0    0.00    0.00     120         std::_Rb_tree_node<std::pair<std::string const, std::string> >* std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_create_node<std::pair<std::string const, std::string> const&>(std::pair<std::string const, std::string> const&) [127]
-----------------------------------------------
                0.00    0.00      80/80          Mixture::load(std::string&) [95]
[128]    0.0    0.00    0.00      80         std::vector<Component, std::allocator<Component> >::push_back(Component const&) [128]
-----------------------------------------------
                0.00    0.00      76/76          Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int) [3]
[129]    0.0    0.00    0.00      76         int minimum<int>(int, int) [129]
-----------------------------------------------
                0.00    0.00      48/48          Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int) [3]
[130]    0.0    0.00    0.00      48         Segment::setInitialDistances(double, double) [130]
-----------------------------------------------
                0.00    0.00       3/27          std::_Rb_tree<std::string, std::pair<std::string const, boost::program_options::variable_value>, std::_Select1st<std::pair<std::string const, boost::program_options::variable_value> >, std::less<std::string>, std::allocator<std::pair<std::string const, boost::program_options::variable_value> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, boost::program_options::variable_value> >*) [184]
                0.00    0.00      24/27          parseCommandLineInput(int, char**) [67]
[131]    0.0    0.00    0.00      27         boost::detail::sp_counted_base::release() [131]
-----------------------------------------------
                0.00    0.00       8/16          parsePDBFile(std::string&) [23]
                0.00    0.00       8/16          std::_Sp_counted_base<(__gnu_cxx::_Lock_policy)2>::_M_release() [22]
[132]    0.0    0.00    0.00      16         lcb::Model::~Model() [132]
-----------------------------------------------
                                 432             std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, std::string> >*) [133]
                0.00    0.00       1/9           parseCommandLineInput(int, char**) [67]
                0.00    0.00       8/9           parsePDBFile(std::string&) [23]
[133]    0.0    0.00    0.00       9+432     std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, std::string> >*) [133]
                                 432             std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, std::string> >*) [133]
-----------------------------------------------
                0.00    0.00       8/8           parsePDBFile(std::string&) [23]
[134]    0.0    0.00    0.00       8         checkFile(std::string&) [134]
-----------------------------------------------
                0.00    0.00       1/8           Protein::initializeCodeLengthMatrices(int) [97]
                0.00    0.00       7/8           loadIdealModels() [24]
[135]    0.0    0.00    0.00       8         IdealModel::~IdealModel() [135]
-----------------------------------------------
                0.00    0.00       8/8           parseCommandLineInput(int, char**) [67]
[136]    0.0    0.00    0.00       8         boost::program_options::typed_value<std::string, char>* boost::program_options::value<std::string>(std::string*) [136]
-----------------------------------------------
                0.00    0.00       8/8           parsePDBFile(std::string&) [23]
[137]    0.0    0.00    0.00       8         lcb::GenericData::getIdentifier() const [137]
-----------------------------------------------
                0.00    0.00       8/8           std::_Sp_counted_base<(__gnu_cxx::_Lock_policy)2>::_M_release() [22]
[138]    0.0    0.00    0.00       8         std::_Sp_counted_ptr_inplace<lcb::Model, std::allocator<lcb::Model>, (__gnu_cxx::_Lock_policy)2>::_M_dispose() [138]
-----------------------------------------------
                0.00    0.00       8/8           parsePDBFile(std::string&) [23]
[139]    0.0    0.00    0.00       8         std::_Sp_counted_ptr_inplace<lcb::Model, std::allocator<lcb::Model>, (__gnu_cxx::_Lock_policy)2>::_M_get_deleter(std::type_info const&) [139]
-----------------------------------------------
                0.00    0.00       8/8           Protein::Protein(lcb::ProteinStructure*, std::string&) [54]
[140]    0.0    0.00    0.00       8         void std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >::_M_emplace_back_aux<lcb::Point<double> const&>(lcb::Point<double> const&) [140]
-----------------------------------------------
                0.00    0.00       8/8           Protein::initializeCodeLengthMatrices(int) [97]
[141]    0.0    0.00    0.00       8         void std::vector<std::vector<OptimalFit, std::allocator<OptimalFit> >, std::allocator<std::vector<OptimalFit, std::allocator<OptimalFit> > > >::_M_emplace_back_aux<std::vector<OptimalFit, std::allocator<OptimalFit> > const&>(std::vector<OptimalFit, std::allocator<OptimalFit> > const&) [141]
                0.00    0.00     608/5852        OptimalFit::OptimalFit(OptimalFit const&) [116]
-----------------------------------------------
                0.00    0.00       8/8           boost::program_options::basic_parsed_options<char> boost::program_options::parse_command_line<char>(int, char const* const*, boost::program_options::options_description const&, int, boost::function1<std::pair<std::string, std::string>, std::string const&>) [84]
[142]    0.0    0.00    0.00       8         void std::vector<std::string, std::allocator<std::string> >::_M_emplace_back_aux<std::string>(std::string&&) [142]
-----------------------------------------------
                0.00    0.00       8/8           parsePDBFile(std::string&) [23]
[143]    0.0    0.00    0.00       8         std::_Rb_tree<std::string, std::pair<std::string const, char>, std::_Select1st<std::pair<std::string const, char> >, std::less<std::string>, std::allocator<std::pair<std::string const, char> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, char> >*) [143]
-----------------------------------------------
                0.00    0.00       8/8           parsePDBFile(std::string&) [23]
[144]    0.0    0.00    0.00       8         std::_Rb_tree<std::string, std::pair<std::string const, int>, std::_Select1st<std::pair<std::string const, int> >, std::less<std::string>, std::allocator<std::pair<std::string const, int> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, int> >*) [144]
-----------------------------------------------
                0.00    0.00       7/7           loadIdealModels() [24]
[145]    0.0    0.00    0.00       7         std::vector<IdealModel, std::allocator<IdealModel> >::push_back(IdealModel const&) [145]
-----------------------------------------------
                0.00    0.00       1/5           Protein::~Protein() [176]
                0.00    0.00       1/5           Protein::Protein(lcb::ProteinStructure*, std::string&) [54]
                0.00    0.00       3/5           boost::program_options::basic_parsed_options<char> boost::program_options::parse_command_line<char>(int, char const* const*, boost::program_options::options_description const&, int, boost::function1<std::pair<std::string, std::string>, std::string const&>) [84]
[146]    0.0    0.00    0.00       5         std::vector<std::string, std::allocator<std::string> >::~vector() [146]
-----------------------------------------------
                0.00    0.00       4/4           void std::vector<IdealModel, std::allocator<IdealModel> >::_M_emplace_back_aux<IdealModel const&>(IdealModel const&) [150]
[147]    0.0    0.00    0.00       4         IdealModel::IdealModel(IdealModel const&) [147]
-----------------------------------------------
                0.00    0.00       4/4           parseCommandLineInput(int, char**) [67]
[148]    0.0    0.00    0.00       4         boost::program_options::typed_value<int, char>* boost::program_options::value<int>(int*) [148]
-----------------------------------------------
                0.00    0.00       1/4           parseCommandLineInput(int, char**) [67]
                0.00    0.00       3/4           boost::program_options::basic_parsed_options<char> boost::program_options::parse_command_line<char>(int, char const* const*, boost::program_options::options_description const&, int, boost::function1<std::pair<std::string, std::string>, std::string const&>) [84]
[149]    0.0    0.00    0.00       4         boost::function1<std::pair<std::string, std::string>, std::string const&>::clear() [149]
-----------------------------------------------
                0.00    0.00       4/4           loadIdealModels() [24]
[150]    0.0    0.00    0.00       4         void std::vector<IdealModel, std::allocator<IdealModel> >::_M_emplace_back_aux<IdealModel const&>(IdealModel const&) [150]
                0.00    0.00       4/4           IdealModel::IdealModel(IdealModel const&) [147]
-----------------------------------------------
                0.00    0.00       1/4           parseCommandLineInput(int, char**) [67]
                0.00    0.00       3/4           boost::program_options::basic_command_line_parser<char>::run() [85]
[151]    0.0    0.00    0.00       4         std::vector<boost::program_options::basic_option<char>, std::allocator<boost::program_options::basic_option<char> > >::~vector() [151]
-----------------------------------------------
                0.00    0.00       4/4           Protein::computeOptimalSegmentation(int) [175]
[152]    0.0    0.00    0.00       4         void std::vector<int, std::allocator<int> >::_M_emplace_back_aux<int const&>(int const&) [152]
-----------------------------------------------
                0.00    0.00       1/2           Protein::translateProteinToOrigin(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [174]
                0.00    0.00       1/2           Protein::Protein(lcb::ProteinStructure*, std::string&) [54]
[153]    0.0    0.00    0.00       2         writeToFile(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&, char const*) [153]
-----------------------------------------------
                0.00    0.00       2/2           parseCommandLineInput(int, char**) [67]
[154]    0.0    0.00    0.00       2         getPDBFilePath(std::string&) [154]
-----------------------------------------------
                0.00    0.00       2/2           Protein::~Protein() [176]
[155]    0.0    0.00    0.00       2         std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::~vector() [155]
-----------------------------------------------
                0.00    0.00       2/2           Protein::computeOptimalSegmentation(int) [175]
[156]    0.0    0.00    0.00       2         void std::vector<int, std::allocator<int> >::_M_emplace_back_aux<int>(int&&) [156]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [681]
[157]    0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN10IdealModelC2Ev [157]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [681]
[158]    0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN10OptimalFitC2Ev [158]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [681]
[159]    0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN10VonMises3DC2Ev [159]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [681]
[160]    0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN16Superpose3DClass23computeRotationalCenterEv [160]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [681]
[161]    0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN6NormalC2Ev [161]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [681]
[162]    0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN7MessageC2Ev [162]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [681]
[163]    0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN7MixtureC2Ev [163]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [681]
[164]    0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN7ProteinC2Ev [164]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [681]
[165]    0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN7SegmentC2EiiRSt6vectorIN3lcb5PointIdEESaIS3_EERS0_ISt5arrayIdLm3EESaIS8_EE [165]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [681]
[166]    0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN9ComponentC2Ev [166]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [681]
[167]    0.0    0.00    0.00       1         _GLOBAL__sub_I_initialize_components_from_file [167]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [681]
[168]    0.0    0.00    0.00       1         _GLOBAL__sub_I_main [168]
-----------------------------------------------
                0.00    0.00       1/1           Protein::translateProteinToOrigin(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [174]
[169]    0.0    0.00    0.00       1         std::ostream& lcb::operator<< <double>(std::ostream&, lcb::Point<double> const&) [169]
-----------------------------------------------
                0.00    0.00       1/1           Protein::printCodeLengthMatrix(int) [173]
[170]    0.0    0.00    0.00       1         char* boost::detail::lcast_put_unsigned<std::char_traits<char>, unsigned int, char>(unsigned int, char*) [170]
-----------------------------------------------
                0.00    0.00       1/1           assignSecondaryStructure(std::string, std::string, int) [1]
[171]    0.0    0.00    0.00       1         Mixture::~Mixture() [171]
-----------------------------------------------
                0.00    0.00       1/1           Protein::Protein(lcb::ProteinStructure*, std::string&) [54]
[172]    0.0    0.00    0.00       1         Protein::checkChainBreak(std::string&, std::vector<lcb::Atom, std::allocator<lcb::Atom> >&) [172]
-----------------------------------------------
                0.00    0.00       1/1           Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int) [3]
[173]    0.0    0.00    0.00       1         Protein::printCodeLengthMatrix(int) [173]
                0.00    0.00       1/1           char* boost::detail::lcast_put_unsigned<std::char_traits<char>, unsigned int, char>(unsigned int, char*) [170]
-----------------------------------------------
                0.00    0.00       1/1           Protein::Protein(lcb::ProteinStructure*, std::string&) [54]
[174]    0.0    0.00    0.00       1         Protein::translateProteinToOrigin(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [174]
                0.00    0.00       1/2           writeToFile(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&, char const*) [153]
                0.00    0.00       1/1           std::ostream& lcb::operator<< <double>(std::ostream&, lcb::Point<double> const&) [169]
-----------------------------------------------
                0.00    0.00       1/1           Protein::compressUsingIdealModels(Mixture&, int) [2]
[175]    0.0    0.00    0.00       1         Protein::computeOptimalSegmentation(int) [175]
                0.00    0.00       4/4           void std::vector<int, std::allocator<int> >::_M_emplace_back_aux<int const&>(int const&) [152]
                0.00    0.00       2/2           void std::vector<int, std::allocator<int> >::_M_emplace_back_aux<int>(int&&) [156]
                0.00    0.00       1/1           std::vector<int, std::allocator<int> >::operator=(std::vector<int, std::allocator<int> > const&) [182]
-----------------------------------------------
                0.00    0.00       1/1           assignSecondaryStructure(std::string, std::string, int) [1]
[176]    0.0    0.00    0.00       1         Protein::~Protein() [176]
                0.00    0.00       2/2           std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::~vector() [155]
                0.00    0.00       1/5           std::vector<std::string, std::allocator<std::string> >::~vector() [146]
-----------------------------------------------
                0.00    0.00       1/1           Protein::initializeCodeLengthMatrices(int) [97]
[177]    0.0    0.00    0.00       1         std::vector<OptimalFit, std::allocator<OptimalFit> >::~vector() [177]
-----------------------------------------------
                0.00    0.00       1/1           void std::vector<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >, std::allocator<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > > >::_M_emplace_back_aux<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > const&>(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > const&) [179]
[178]    0.0    0.00    0.00       1         std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >::vector(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > const&) [178]
-----------------------------------------------
                0.00    0.00       1/1           Protein::Protein(lcb::ProteinStructure*, std::string&) [54]
[179]    0.0    0.00    0.00       1         void std::vector<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >, std::allocator<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > > >::_M_emplace_back_aux<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > const&>(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > const&) [179]
                0.00    0.00       1/1           std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >::vector(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > const&) [178]
-----------------------------------------------
                0.00    0.00       1/1           Protein::computeSphericalTransformation() [81]
[180]    0.0    0.00    0.00       1         void std::vector<std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > >, std::allocator<std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > > > >::_M_emplace_back_aux<std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > > const&>(std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > > const&) [180]
-----------------------------------------------
                0.00    0.00       1/1           Protein::Protein(lcb::ProteinStructure*, std::string&) [54]
[181]    0.0    0.00    0.00       1         void std::vector<std::string, std::allocator<std::string> >::_M_emplace_back_aux<std::string const&>(std::string const&) [181]
-----------------------------------------------
                0.00    0.00       1/1           Protein::computeOptimalSegmentation(int) [175]
[182]    0.0    0.00    0.00       1         std::vector<int, std::allocator<int> >::operator=(std::vector<int, std::allocator<int> > const&) [182]
-----------------------------------------------
                0.00    0.00       1/1           parseCommandLineInput(int, char**) [67]
[183]    0.0    0.00    0.00       1         std::_Rb_tree<std::string, std::string, std::_Identity<std::string>, std::less<std::string>, std::allocator<std::string> >::_M_erase(std::_Rb_tree_node<std::string>*) [183]
-----------------------------------------------
                                   3             std::_Rb_tree<std::string, std::pair<std::string const, boost::program_options::variable_value>, std::_Select1st<std::pair<std::string const, boost::program_options::variable_value> >, std::less<std::string>, std::allocator<std::pair<std::string const, boost::program_options::variable_value> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, boost::program_options::variable_value> >*) [184]
                0.00    0.00       1/1           parseCommandLineInput(int, char**) [67]
[184]    0.0    0.00    0.00       1+3       std::_Rb_tree<std::string, std::pair<std::string const, boost::program_options::variable_value>, std::_Select1st<std::pair<std::string const, boost::program_options::variable_value> >, std::less<std::string>, std::allocator<std::pair<std::string const, boost::program_options::variable_value> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, boost::program_options::variable_value> >*) [184]
                0.00    0.00       3/27          boost::detail::sp_counted_base::release() [131]
                                   3             std::_Rb_tree<std::string, std::pair<std::string const, boost::program_options::variable_value>, std::_Select1st<std::pair<std::string const, boost::program_options::variable_value> >, std::less<std::string>, std::allocator<std::pair<std::string const, boost::program_options::variable_value> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, boost::program_options::variable_value> >*) [184]
-----------------------------------------------
                0.00    0.00       1/1           parseCommandLineInput(int, char**) [67]
[185]    0.0    0.00    0.00       1         std::basic_ostream<char, std::char_traits<char> >& std::operator<< <std::char_traits<char> >(std::basic_ostream<char, std::char_traits<char> >&, char const*) [clone .constprop.477] [185]
-----------------------------------------------
                                 112             std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_copy(std::_Rb_tree_node<std::pair<std::string const, std::string> > const*, std::_Rb_tree_node<std::pair<std::string const, std::string> >*) [580]
[580]    0.0    0.00    0.00       0+112     std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_copy(std::_Rb_tree_node<std::pair<std::string const, std::string> > const*, std::_Rb_tree_node<std::pair<std::string const, std::string> >*) [580]
                0.00    0.00     120/120         std::_Rb_tree_node<std::pair<std::string const, std::string> >* std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_create_node<std::pair<std::string const, std::string> const&>(std::pair<std::string const, std::string> const&) [127]
                                 112             std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_copy(std::_Rb_tree_node<std::pair<std::string const, std::string> > const*, std::_Rb_tree_node<std::pair<std::string const, std::string> >*) [580]
-----------------------------------------------

 This table describes the call tree of the program, and was sorted by
 the total amount of time spent in each function and its children.

 Each entry in this table consists of several lines.  The line with the
 index number at the left hand margin lists the current function.
 The lines above it list the functions that called this function,
 and the lines below it list the functions this one called.
 This line lists:
     index	A unique number given to each element of the table.
		Index numbers are sorted numerically.
		The index number is printed next to every function name so
		it is easier to look up where the function is in the table.

     % time	This is the percentage of the `total' time that was spent
		in this function and its children.  Note that due to
		different viewpoints, functions excluded by options, etc,
		these numbers will NOT add up to 100%.

     self	This is the total amount of time spent in this function.

     children	This is the total amount of time propagated into this
		function by its children.

     called	This is the number of times the function was called.
		If the function called itself recursively, the number
		only includes non-recursive calls, and is followed by
		a `+' and the number of recursive calls.

     name	The name of the current function.  The index number is
		printed after it.  If the function is a member of a
		cycle, the cycle number is printed between the
		function's name and the index number.


 For the function's parents, the fields have the following meanings:

     self	This is the amount of time that was propagated directly
		from the function into this parent.

     children	This is the amount of time that was propagated from
		the function's children into this parent.

     called	This is the number of times this parent called the
		function `/' the total number of times the function
		was called.  Recursive calls to the function are not
		included in the number after the `/'.

     name	This is the name of the parent.  The parent's index
		number is printed after it.  If the parent is a
		member of a cycle, the cycle number is printed between
		the name and the index number.

 If the parents of the function cannot be determined, the word
 `<spontaneous>' is printed in the `name' field, and all the other
 fields are blank.

 For the function's children, the fields have the following meanings:

     self	This is the amount of time that was propagated directly
		from the child into the function.

     children	This is the amount of time that was propagated from the
		child's children to the function.

     called	This is the number of times the function called
		this child `/' the total number of times the child
		was called.  Recursive calls by the child are not
		listed in the number after the `/'.

     name	This is the name of the child.  The child's index
		number is printed after it.  If the child is a
		member of a cycle, the cycle number is printed
		between the name and the index number.

 If there are any cycles (circles) in the call graph, there is an
 entry for the cycle-as-a-whole.  This entry shows who called the
 cycle (as parents) and the members of the cycle (as children.)
 The `+' recursive calls entry shows the number of function calls that
 were internal to the cycle, and the calls entry for each member shows,
 for that member, how many times it was called from other members of
 the cycle.

Copyright (C) 2012 Free Software Foundation, Inc.

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.

Index by function name

 [157] _GLOBAL__sub_I__ZN10IdealModelC2Ev (IdealModel.cpp) [82] lcb::Matrix<double>::Matrix(int, int) [76] Component::computeSecondDerivative(double)
 [158] _GLOBAL__sub_I__ZN10OptimalFitC2Ev (OptimalFit.cpp) [108] lcb::Matrix<double>::Matrix() [6] Component::conflate(Component&)
 [159] _GLOBAL__sub_I__ZN10VonMises3DC2Ev (VonMises3D.cpp) [58] lcb::Matrix<double>::~Matrix() [90] Component::generate(int)
 [160] _GLOBAL__sub_I__ZN16Superpose3DClass23computeRotationalCenterEv (Superpose3D.cpp) [32] lcb::Matrix<double>::operator=(lcb::Matrix<double> const&) [66] Component::Component(std::array<double, 2ul>&, double, int)
 [161] _GLOBAL__sub_I__ZN6NormalC2Ev (Normal.cpp) [29] lcb::Matrix<double>::operator*(lcb::Matrix<double> const&) [69] Component::operator=(Component const&)
 [162] _GLOBAL__sub_I__ZN7MessageC2Ev (Message.cpp) [42] lcb::Vector<double>::angleBetween(lcb::Vector<double> const&, lcb::Vector<double> const&) [120] OptimalFit::getMessageLength() const
 [163] _GLOBAL__sub_I__ZN7MixtureC2Ev (Mixture.cpp) [57] lcb::Vector<double>::crossProduct(lcb::Vector<double> const&, lcb::Vector<double> const&) [137] lcb::GenericData::getIdentifier() const
 [164] _GLOBAL__sub_I__ZN7ProteinC2Ev (Protein.cpp) [9] lcb::Vector<double>::normalize() [26] lcb::Matrix<double>::operator*(lcb::Vector<double> const&) const
 [165] _GLOBAL__sub_I__ZN7SegmentC2EiiRSt6vectorIN3lcb5PointIdEESaIS3_EERS0_ISt5arrayIdLm3EESaIS8_EE (Segment.cpp) [35] lcb::Matrix<double> lcb::geometry::rotationMatrix<double>(lcb::Vector<double> const&, double) (Geometry.h) [44] lcb::Vector<double>::normalize_copy() const
 [166] _GLOBAL__sub_I__ZN9ComponentC2Ev (Component.cpp) [25] lcb::Point<double> lcb::geometry::transform<double>(lcb::Point<double> const&, lcb::Matrix<double> const&) (Geometry.h) [11] lcb::Vector<double>::l2Norm() const
 [167] _GLOBAL__sub_I_initialize_components_from_file (Support.cpp) [38] lcb::Point<double> lcb::geometry::transform<double>(lcb::Point<double> const&, lcb::Matrix<double> const&) (Geometry.h) [59] lcb::Vector<double>::operator*(lcb::Vector<double> const&) const
 [168] _GLOBAL__sub_I_main (main.cpp) [169] std::ostream& lcb::operator<< <double>(std::ostream&, lcb::Point<double> const&) [73] std::_Rb_tree<std::string, std::pair<std::string const, boost::program_options::variable_value>, std::_Select1st<std::pair<std::string const, boost::program_options::variable_value> >, std::less<std::string>, std::allocator<std::pair<std::string const, boost::program_options::variable_value> > >::find(std::string const&) const
  [63] extractName(std::string&) [125] boost::char_separator<char, std::char_traits<char> >::char_separator(boost::char_separator<char, std::char_traits<char> > const&) [22] std::_Sp_counted_base<(__gnu_cxx::_Lock_policy)2>::_M_release()
 [153] writeToFile(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&, char const*) [123] boost::char_separator<char, std::char_traits<char> >::~char_separator() [77] std::_Sp_counted_base<(__gnu_cxx::_Lock_policy)2>::~_Sp_counted_base()
  [23] parsePDBFile(std::string&) [122] bool boost::char_separator<char, std::char_traits<char> >::operator()<__gnu_cxx::__normal_iterator<char const*, std::string>, std::string>(__gnu_cxx::__normal_iterator<char const*, std::string>&, __gnu_cxx::__normal_iterator<char const*, std::string>, std::string&) [138] std::_Sp_counted_ptr_inplace<lcb::Model, std::allocator<lcb::Model>, (__gnu_cxx::_Lock_policy)2>::_M_dispose()
  [50] point2vector(lcb::Point<double>&) [124] boost::token_iterator<boost::char_separator<char, std::char_traits<char> >, __gnu_cxx::__normal_iterator<char const*, std::string>, std::string>::~token_iterator() [139] std::_Sp_counted_ptr_inplace<lcb::Model, std::allocator<lcb::Model>, (__gnu_cxx::_Lock_policy)2>::_M_get_deleter(std::type_info const&)
 [154] getPDBFilePath(std::string&) [84] boost::program_options::basic_parsed_options<char> boost::program_options::parse_command_line<char>(int, char const* const*, boost::program_options::options_description const&, int, boost::function1<std::pair<std::string, std::string>, std::string const&>) [70] void std::__uninitialized_default_n_1<false>::__uninit_default_n<lcb::Vector<double>*, unsigned long>(lcb::Vector<double>*, unsigned long)
  [24] loadIdealModels()      [88] boost::program_options::error_with_option_name::set_option_name(std::string const&) [150] void std::vector<IdealModel, std::allocator<IdealModel> >::_M_emplace_back_aux<IdealModel const&>(IdealModel const&)
   [8] convertToCartesian(double, double, double) [85] boost::program_options::basic_command_line_parser<char>::run() [145] std::vector<IdealModel, std::allocator<IdealModel> >::push_back(IdealModel const&)
   [7] convertToSpherical(lcb::Point<double>&) [136] boost::program_options::typed_value<std::string, char>* boost::program_options::value<std::string>(std::string*) [74] std::vector<IdealModel, std::allocator<IdealModel> >::~vector()
  [86] generateRandomWeights(int, double) [148] boost::program_options::typed_value<int, char>* boost::program_options::value<int>(int*) [177] std::vector<OptimalFit, std::allocator<OptimalFit> >::~vector()
  [10] convertToCanonicalForm(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [131] boost::detail::sp_counted_base::release() [17] void std::vector<Component, std::allocator<Component> >::_M_emplace_back_aux<Component const&>(Component const&)
   [1] assignSecondaryStructure(std::string, std::string, int) [170] char* boost::detail::lcast_put_unsigned<std::char_traits<char>, unsigned int, char>(unsigned int, char*) [128] std::vector<Component, std::allocator<Component> >::push_back(Component const&)
 [129] int minimum<int>(int, int) [149] boost::function1<std::pair<std::string, std::string>, std::string const&>::clear() [55] std::vector<lcb::Atom, std::allocator<lcb::Atom> >::~vector()
 [134] checkFile(std::string&) [64] Normal::negativeLogLikelihood(std::vector<double, std::allocator<double> >&) [140] void std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >::_M_emplace_back_aux<lcb::Point<double> const&>(lcb::Point<double> const&)
  [80] IdealModel::getResidues(int) [106] Normal::density(double) [178] std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >::vector(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > const&)
 [112] IdealModel::getStructure() [109] Normal::Normal(double, double) [91] std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >::~vector()
 [113] IdealModel::getName() [107] Message::encodeUsingNormalModel(double, Normal&) [52] std::vector<lcb::Vector<double>, std::allocator<lcb::Vector<double> > >::_M_default_append(unsigned long)
 [114] IdealModel::setLength(int) [121] Message::encodeUsingSphereModel(double, Normal&) [46] std::vector<lcb::Vector<double>, std::allocator<lcb::Vector<double> > >::~vector()
 [117] IdealModel::IdealModel() [79] Message::encodeUsingLogStarModel(double) [151] std::vector<boost::program_options::basic_option<char>, std::allocator<boost::program_options::basic_option<char> > >::~vector()
  [34] IdealModel::IdealModel(lcb::ProteinStructure*, std::string) [14] Message::encodeUsingMixtureModel(std::array<double, 2ul>&, Mixture&) [141] void std::vector<std::vector<OptimalFit, std::allocator<OptimalFit> >, std::allocator<std::vector<OptimalFit, std::allocator<OptimalFit> > > >::_M_emplace_back_aux<std::vector<OptimalFit, std::allocator<OptimalFit> > const&>(std::vector<OptimalFit, std::allocator<OptimalFit> > const&)
 [147] IdealModel::IdealModel(IdealModel const&) [110] Message::Message() [78] std::vector<std::vector<OptimalFit, std::allocator<OptimalFit> >, std::allocator<std::vector<OptimalFit, std::allocator<OptimalFit> > > >::~vector()
 [119] IdealModel::IdealModel(int, std::string) [89] Mixture::initialize2() [179] void std::vector<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >, std::allocator<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > > >::_M_emplace_back_aux<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > const&>(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > > const&)
 [135] IdealModel::~IdealModel() [15] Mixture::probability(std::array<double, 2ul>&) [180] void std::vector<std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > >, std::allocator<std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > > > >::_M_emplace_back_aux<std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > > const&>(std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > > const&)
  [65] IdealModel::operator=(IdealModel const&) [95] Mixture::load(std::string&) [40] void std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::_M_emplace_back_aux<std::vector<double, std::allocator<double> > const&>(std::vector<double, std::allocator<double> > const&)
  [62] OptimalFit::OptimalFit(IdealModel&, double) [5] Mixture::conflate(Component&) [155] std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::~vector()
 [116] OptimalFit::OptimalFit(OptimalFit const&) [21] Mixture::Mixture(int, std::vector<Component, std::allocator<Component> >&, std::vector<double, std::allocator<double> >&) [28] std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >::operator=(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > > const&)
 [118] OptimalFit::OptimalFit() [111] Mixture::Mixture() [181] void std::vector<std::string, std::allocator<std::string> >::_M_emplace_back_aux<std::string const&>(std::string const&)
  [49] OptimalFit::operator=(OptimalFit const&) [171] Mixture::~Mixture() [142] void std::vector<std::string, std::allocator<std::string> >::_M_emplace_back_aux<std::string>(std::string&&)
 [115] OptimalFit::operator<(OptimalFit const&) [172] Protein::checkChainBreak(std::string&, std::vector<lcb::Atom, std::allocator<lcb::Atom> >&) [146] std::vector<std::string, std::allocator<std::string> >::~vector()
  [19] VonMises3D::density(double, double) [93] Protein::computeTransformation(int, int) [83] void std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > >::_M_emplace_back_aux<std::array<double, 3ul> const&>(std::array<double, 3ul> const&)
  [18] VonMises3D::VonMises3D(std::array<double, 2ul>&, double) [173] Protein::printCodeLengthMatrix(int) [41] std::vector<double, std::allocator<double> >::_M_fill_insert(__gnu_cxx::__normal_iterator<double*, std::vector<double, std::allocator<double> > >, unsigned long, double const&)
  [20] VonMises3D::VonMises3D() [3] Protein::computeCodeLengthMatrix(std::vector<IdealModel, std::allocator<IdealModel> >&, Mixture&, int, int) [56] std::vector<double, std::allocator<double> >::_M_default_append(unsigned long)
  [39] VonMises3D::operator=(VonMises3D const&) [2] Protein::compressUsingIdealModels(Mixture&, int) [51] void std::vector<double, std::allocator<double> >::_M_emplace_back_aux<double const&>(double const&)
  [16] Superpose3DClass::diagonalize() [174] Protein::translateProteinToOrigin(std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&) [92] std::vector<double, std::allocator<double> >::resize(unsigned long)
  [48] Superpose3DClass::transformVectors(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&) [175] Protein::computeOptimalSegmentation(int) [126] std::vector<double, std::allocator<double> >::push_back(double const&)
  [61] Superpose3DClass::updateCenterOfMasses() [98] Protein::computeSuccessiveDistances() [60] std::vector<double, std::allocator<double> >::vector(unsigned long, double const&, std::allocator<double> const&)
  [72] Superpose3DClass::computeQuaternionMatrix() [97] Protein::initializeCodeLengthMatrices(int) [13] std::vector<double, std::allocator<double> >::vector(std::vector<double, std::allocator<double> > const&)
  [71] Superpose3DClass::getSufficientStatistics() [81] Protein::computeSphericalTransformation() [31] std::vector<double, std::allocator<double> >::operator=(std::vector<double, std::allocator<double> > const&)
  [75] Superpose3DClass::getRMSD() [94] Protein::computeMessageLengthUsingNullModel(Mixture&) [152] void std::vector<int, std::allocator<int> >::_M_emplace_back_aux<int const&>(int const&)
  [45] Superpose3DClass::Superpose3DClass(std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&, std::vector<std::vector<double, std::allocator<double> >, std::allocator<std::vector<double, std::allocator<double> > > >&) [96] Protein::computeMessageLengthUsingSphereModel() [156] void std::vector<int, std::allocator<int> >::_M_emplace_back_aux<int>(int&&)
  [36] Superpose3DClass::Superpose3DClass(std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&) [54] Protein::Protein(lcb::ProteinStructure*, std::string&) [182] std::vector<int, std::allocator<int> >::operator=(std::vector<int, std::allocator<int> > const&)
  [12] Superpose3DClass::Superpose3DClass(suffStatClass&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&) [176] Protein::~Protein() [183] std::_Rb_tree<std::string, std::string, std::_Identity<std::string>, std::less<std::string>, std::allocator<std::string> >::_M_erase(std::_Rb_tree_node<std::string>*)
  [30] lcb::GenericData::~GenericData() [33] Segment::fitNullModel(Mixture&) [184] std::_Rb_tree<std::string, std::pair<std::string const, boost::program_options::variable_value>, std::_Select1st<std::pair<std::string const, boost::program_options::variable_value> >, std::less<std::string>, std::allocator<std::pair<std::string const, boost::program_options::variable_value> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, boost::program_options::variable_value> >*)
  [53] lcb::Chain::~Chain()    [4] Segment::fitIdealModel(IdealModel&, Mixture&, int) [127] std::_Rb_tree_node<std::pair<std::string const, std::string> >* std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_create_node<std::pair<std::string const, std::string> const&>(std::pair<std::string const, std::string> const&)
  [87] lcb::Model::~Model()  [130] Segment::setInitialDistances(double, double) [133] std::_Rb_tree<std::string, std::pair<std::string const, std::string>, std::_Select1st<std::pair<std::string const, std::string> >, std::less<std::string>, std::allocator<std::pair<std::string const, std::string> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, std::string> >*)
 [132] lcb::Model::~Model()   [27] Segment::getCurrentMeanAndDirection(std::pair<std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >, lcb::Matrix<double> >&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&, int) [143] std::_Rb_tree<std::string, std::pair<std::string const, char>, std::_Select1st<std::pair<std::string const, char> >, std::less<std::string>, std::allocator<std::pair<std::string const, char> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, char> >*)
  [37] lcb::Matrix<double>::Matrix(lcb::Matrix<double> const&) [68] Segment::Segment(int, int, std::vector<lcb::Point<double>, std::allocator<lcb::Point<double> > >&, std::vector<std::array<double, 3ul>, std::allocator<std::array<double, 3ul> > >&) [144] std::_Rb_tree<std::string, std::pair<std::string const, int>, std::_Select1st<std::pair<std::string const, int> >, std::less<std::string>, std::allocator<std::pair<std::string const, int> > >::_M_erase(std::_Rb_tree_node<std::pair<std::string const, int> >*)
  [43] lcb::Matrix<double>::Matrix(int) [47] Component::likelihood(std::array<double, 2ul>&) [185] std::basic_ostream<char, std::char_traits<char> >& std::operator<< <std::char_traits<char> >(std::basic_ostream<char, std::char_traits<char> >&, char const*) [clone .constprop.477] (ostream)
