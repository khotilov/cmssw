#ifndef ECALDETID_ECALCONTAINER_H
#define ECALDETID_ECALCONTAINER_H

#include "DataFormats/DetId/interface/DetId.h"
#include <vector>
#include <utility>
#include <algorithm>


/* a generic container for ecal items
 * provides access by hashedIndex and by DetId...
 */

template<typename DetId, typename T>
class EcalContainer {

        public:

                typedef EcalContainer<DetId, T> self;
                typedef T Item;
                typedef Item value_type;
                typedef typename std::vector<Item> Items; 
                typedef typename std::vector<Item>::const_iterator const_iterator;
                typedef typename std::vector<Item>::iterator iterator;

                void insert(std::pair<uint32_t, Item> const &a) {
                        (*this)[a.first] = a.second;
                }

                inline const Item & item(size_t hashid) const {
                        return m_items[hashid];
                }

                inline const Items & items() const {
                        return m_items;
                }

                inline Item & operator[](::DetId id) {
                        if (m_items.empty()) {
                                m_items.resize(DetId::SIZE_HASH);
                        }
                        static Item dummy;
                        DetId ib(id.rawId());
                        if ( !isValidId(ib) ) return dummy;
                        return m_items[ib.hashedIndex()];
                }

                inline Item const & operator[](::DetId id) const {
                        static Item dummy;
                        DetId ib(id.rawId());
                        if ( !isValidId(ib) ) return dummy;
                        return m_items[ib.hashedIndex()];
                }

                inline Item & operator[](uint32_t rawId) {
                        if (m_items.empty()) {
                                m_items.resize(DetId::SIZE_HASH);
                        }
                        static Item dummy;
                        DetId id(rawId);
                        if ( !isValidId(id) ) return dummy;
                        return m_items[id.hashedIndex()];
                }

                inline Item const & operator[](uint32_t rawId) const {
                        static Item dummy;
                        DetId id(rawId);
                        if ( !isValidId(id) ) return dummy;
                        return m_items[id.hashedIndex()];
                }

                inline const_iterator find(uint32_t rawId) const {
                        DetId ib(rawId);
                        if ( !isValidId(ib) ) return m_items.end();
                        return m_items.begin() + ib.hashedIndex();
                }

                inline const_iterator begin() const {
                        return m_items.begin();
                }

                inline const_iterator end() const {
                        return m_items.end();
                }

        private:

                // not protected on EB <--> EE swap -- FIXME?
                inline bool isValidId(const DetId id) const {
                        return id.det() == ::DetId::Ecal;
                };

                std::vector<Item> m_items;

};



#endif // ECALCONTAINER
