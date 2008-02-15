#ifndef GeneratorInterface_LHEInterface_LHEReader_h
#define GeneratorInterface_LHEInterface_LHEReader_h

#include <string>
#include <vector>
#include <memory>

#include <boost/shared_ptr.hpp>

#include "FWCore/ParameterSet/interface/ParameterSet.h"

namespace lhef {

class XMLDocument;
class LHECommon;
class LHEEvent;

class LHEReader {
    public:
	LHEReader(const edm::ParameterSet &params);
	~LHEReader();

	boost::shared_ptr<LHEEvent> next();

    private:
	class Source;
	class FileSource;
	class XMLHandler;

	const std::vector<std::string>	fileURLs;
	unsigned int			firstEvent;
	unsigned int			curIndex;

	std::auto_ptr<Source>		curSource;
	std::auto_ptr<XMLDocument>	curDoc;
	boost::shared_ptr<LHECommon>	curCommon;
	std::auto_ptr<XMLHandler>	handler;
};

} // namespace lhef

#endif // GeneratorInterface_LHEInterface_LHEReader_h
