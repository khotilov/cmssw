#ifndef DDLSAX2HANDLER_H
#define DDLSAX2HANDLER_H

// ---------------------------------------------------------------------------
//  Includes
// ---------------------------------------------------------------------------
#include "DetectorDescription/Parser/interface/DDLSAX2Handler.h"

// Xerces C++ Dependencies
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/sax2/DefaultHandler.hpp>

#include <string>
#include <vector>
#include <iostream>

namespace std{} using namespace std;
namespace xercesc_2_3{} using namespace xercesc_2_3;

///  DDLSAX2Handler inherits from Xerces C++ DefaultHandler.
/** @class DDLSAX2Handler
 * @author Michael Case
 *
 *  DDLSAX2Handler.h  -  description
 *  -------------------
 *  begin: Mon Oct 22 2001
 *  email: case@ucdhep.ucdavis.edu
 *
 *  The DefaultHandler of Xerces C++ provides an interface to the SAX2 event
 *  driven processing of XML documents.  It does so by providing methods which
 *  are redefined by the inheriting class (DDLSAX2Handler in this case) to
 *  provide the desired processing for each event.
 *
 *  In this case, we accumulate some statistics.  This class does nothing with
 *  startElement and endElement events.
 *
 */
class DDLSAX2Handler : public DefaultHandler
{

 public:

  // -----------------------------------------------------------------------
  //  Constructor and Destructor
  // -----------------------------------------------------------------------

  DDLSAX2Handler();
  ~DDLSAX2Handler();

  // -----------------------------------------------------------------------
  //  Getter methods
  // -----------------------------------------------------------------------
  /// Get the count of elements processed so far.
  unsigned int getElementCount() const
    {
      return elementCount_;
    }
  /// Get the count of attributes processed so far.
  unsigned int getAttrCount() const
    {
      return attrCount_;
    }
  /// Get the count of characters processed so far.
  unsigned int getCharacterCount() const
    {
      return characterCount_;
    }
  /// Did the XML parser see any errors?
  bool getSawErrors() const
    {
      return sawErrors_;
    }
  /// Get the count of spaces processed so far.
  unsigned int getSpaceCount() const
    {
      return spaceCount_;
    }

  // -----------------------------------------------------------------------
  //  Handlers for the SAX ContentHandler interface
  // -----------------------------------------------------------------------

  virtual void startElement(const XMLCh* const uri, const XMLCh* const localname
		    , const XMLCh* const qname, const Attributes& attrs);
  virtual void endElement(const XMLCh* const uri, const XMLCh* const localname
		    , const XMLCh* const qname);
  virtual void characters(const XMLCh* const chars, const unsigned int length);
  virtual void comment (const XMLCh *const chars, const unsigned int length );
  virtual void ignorableWhitespace(const XMLCh* const chars, const unsigned int length);
  virtual void resetDocument();

  // -----------------------------------------------------------------------
  //  Handlers for the SAX ErrorHandler interface
  // -----------------------------------------------------------------------
  virtual void warning(const SAXParseException& exception);
  virtual void error(const SAXParseException& exception);
  virtual void fatalError(const SAXParseException& exception);
  virtual void dumpStats(const string& fname);

 protected:

  // -----------------------------------------------------------------------
  //  Protected data members
  //
  //  attrCount_
  //  characterCount_
  //  elementCount_
  //  spaceCount_
  //      These are just counters that are run upwards based on the input
  //      from the document handlers.
  //
  //  sawErrors
  //      This is set by the error handlers, and is queryable later to
  //      see if any errors occurred.
  // -----------------------------------------------------------------------
  unsigned int    attrCount_;
  unsigned int    characterCount_;
  unsigned int    elementCount_;
  unsigned int    spaceCount_;
  bool            sawErrors_;

  std::string getnmspace(const std::string& fname);
};

#endif
