// $Id: InitMsgCollection.h,v 1.8 2009/09/16 10:44:29 mommsen Exp $
/// @file: InitMsgCollection.h 

#ifndef StorageManager_InitMsgCollection_h
#define StorageManager_InitMsgCollection_h

#include "EventFilter/StorageManager/interface/ConsumerID.h"

#include "IOPool/Streamer/interface/InitMessage.h"

#include "boost/shared_ptr.hpp"
#include "boost/thread/thread.hpp"
#include <vector>
#include <map>
#include <string>
#include <utility>

namespace stor
{

  /**
     This class is used to manage the unique set of INIT messages
     that have been received by the storage manager and will be sent
     to event consumers and written to output streams.

     $Author: mommsen $
     $Revision: 1.8 $
     $Date: 2009/09/16 10:44:29 $
  */

  typedef std::vector<unsigned char> InitMsgBuffer;
  typedef boost::shared_ptr<InitMsgBuffer> InitMsgSharedPtr;

  class InitMsgCollection
  {

  public:

    /**
     * InitMsgCollection constructor.
     */
    InitMsgCollection();

    /**
     * InitMsgCollection destructor.
     */
    ~InitMsgCollection();

    /**
     * Adds the specified INIT message to the collection if it has a unique
     * HLT output module label.
     *
     * If we already have an INIT message with the same output module label
     * as the input INIT message, the duplicate
     * message is *not* added to the collection, and this method returns false.
     *
     * If the output module label inside the INIT message is empty, an
     * exception is thrown.
     *
     * @param initMsgView The INIT message to be added to the collection.
     * @return true if the message was added to the collection, false otherwise.
     * @throws cms::Exception if one of the consistency checks fails.
     */
    bool addIfUnique(InitMsgView const& initMsgView);

    /**
     * Fetches the single INIT message that matches the requested HLT output
     * module label.  If no messages match the request, an empty pointer
     * is returned.
     *
     * If the requested HLT output module label is empty, and there is only
     * one INIT message in the collection, that INIT message is returned.
     * However, if there is more than one INIT message in the collection, and
     * an empty request is passed into this method, an exception will be thrown.
     * (How can we decide which is the best match when we have nothing to work
     * with?)
     *
     * @param requestedOMLabel The HLT output module label of the INIT
     *        message to be returned.
     * @return a pointer to the INIT message that matches.  If no
     *         matching INIT message is found, and empty pointer is returned.
     * @throws cms::Exception if the input HLT output module label string is
     *         empty and there is more than one INIT message in the collection.
     */
    InitMsgSharedPtr getElementForOutputModule(const std::string& requestedOMLabel) const;

    /**
     * Returns a shared pointer to the last element in the collection
     * or an empty pointer if the collection has no elements.
     *
     * @return the last InitMsgSharedPtr in the collection.
     */
    InitMsgSharedPtr getLastElement() const;

    /**
     * Returns a shared pointer to the requested element in the collection
     * or an empty pointer if the requested index if out of bounds.
     *
     * @param index The index of the requested element.
     * @return the InitMsgSharedPtr at the requested index in the collection.
     */
    InitMsgSharedPtr getElementAt(const unsigned int index) const;
 
    /**
     * Returns a shared pointer to all INIT messages in the collection
     * or an empty pointer if no collections are stored.
     *
     * @return the InitMsgSharedPtr at the beginning of the full collection.
     */
    InitMsgSharedPtr getFullCollection() const { return serializedFullSet_; }

    /**
     * Removes all entries from the collection.
     */
    void clear();

    /**
     * Returns the number of unique INIT messages in the collection.
     *
     * @return the integer number of messages.
     */
    int size() const;

    /**
     * Returns the number of identical INIT messages received for the
     * given module name
     *
     * @return the integer number of received INIT messages
     */
    uint32 initMsgCount(const std::string& outputModuleLabel) const;

    /**
     * Returns the maximum number of identical INIT messages received
     * for any output module
     *
     * @return the integer number of maximum received INIT messages
     */
    uint32 maxMsgCount() const;

    /**
     * Returns a string with information on which selections are available.
     *
     * @return the help string.
     */
    std::string getSelectionHelpString() const;

    /**
     * Returns the name of the output module with the specified module ID,
     * or an empty string of the specified module ID is not known.
     *
     * @return the output module label or an empty string
     */
    std::string getOutputModuleName(const uint32 outputModuleId) const;

    /**
     * Creates a single text string from the elements in the specified
     * list of strings.  The specified maximum number of elements are
     * included, however a zero value for the maximum number will include
     * all elements.
     *
     * @param list the list of strings to include (vector of strings);
     * @param maxCount the maximum number of list elements to include.
     * @return the text string with the formatted list elements.
     */
    static std::string stringsToText(Strings const& list,
                                     unsigned int maxCount = 0);

  private:

    /**
     * Adds the specified INIT message to the collection (unconditionally).
     *
     * @param initMsgView The INIT message to add to the collection.
     */
    void add(InitMsgView const& initMsgView);

    typedef std::pair<InitMsgSharedPtr, uint32> InitMsgPtrAndCount;
    typedef std::vector<InitMsgPtrAndCount> InitMsgList;
    InitMsgList initMsgList_;
    InitMsgSharedPtr serializedFullSet_;

    typedef std::map<uint32, std::string> OutModTable;
    OutModTable outModNameTable_;
    mutable boost::mutex listLock_;
  };
}

#endif // StorageManager_InitMsgCollection_h


/// emacs configuration
/// Local Variables: -
/// mode: c++ -
/// c-basic-offset: 2 -
/// indent-tabs-mode: nil -
/// End: -
