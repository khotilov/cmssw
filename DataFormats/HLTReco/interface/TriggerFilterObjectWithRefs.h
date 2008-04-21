#ifndef HLTReco_TriggerFilterObjectWithRefs_h
#define HLTReco_TriggerFilterObjectWithRefs_h

/** \class trigger::TriggerFilterObjectWithRefs
 *
 *  If HLT cuts of intermediate or final HLT filters are satisfied,
 *  instances of this class hold the combination of reconstructed
 *  physics objects (e/gamma/mu/jet/MMet...) satisfying the cuts.
 *
 *  This implementation is not completely space-efficient as some
 *  physics object containers may stay empty. However, the big
 *  advantage is that the solution is generic, i.e., works for all
 *  possible HLT filters. Hence we accept the reasonably small
 *  overhead of empty containers.
 *
 *  $Date: 2007/12/05 14:24:02 $
 *  $Revision: 1.9 $
 *
 *  \author Martin Grunewald
 *
 */

#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include <string>
#include <vector>

namespace trigger
{

  /// Transient book-keeping EDProduct filled by HLTFilter module to
  /// record physics objects firing the filter (never persistet in
  /// production; same functionality but different implementation
  /// compared to the old HLT data model's HLTFilterObjectWithRefs
  /// class)
  class TriggerFilterObjectWithRefs : public TriggerRefsCollections {

  /// data members
  private:
    int path_;
    int module_;
    std::vector<std::string> collectionTags_;

  /// methods
  public:
    /// constructors
    TriggerFilterObjectWithRefs():
      TriggerRefsCollections(),
      path_(-9),
      module_(-9),
      collectionTags_() { }
    TriggerFilterObjectWithRefs(int path, int module):
      TriggerRefsCollections(),
      path_(path),
      module_(module),
      collectionTags_() { }
    /// accessors
    int path() const {return path_;}
    int module() const {return module_;}
    /// tags
    const std::vector<edm::InputTag> getCollectionTags() const {
      const size_type n(collectionTags_.size());
      std::vector<edm::InputTag> tags(n);
      for (size_type i=0; i!=n; ++i) {
	tags[i]=edm::InputTag(collectionTags_[i]);
      }
      return tags;
    }
    void addCollectionTag(const edm::InputTag& collectionTag){
      collectionTags_.push_back(collectionTag.encode());
      return;
    }
  };

}

#endif
