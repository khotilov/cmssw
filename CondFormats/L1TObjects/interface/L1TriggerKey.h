#ifndef CondFormats_L1TObjects_L1TriggerKey_h
#define CondFormats_L1TObjects_L1TriggerKey_h

#include <string>
#include <map>

/* L1 key used to load all other configuration data from offline db.
 * This class is just a proxy to the real data. It will contain mapping from data and record
 * pair to the payload token that could be used to read data. So the use case could be as follows:
 *   1. User read L1TriggerKey for given Tag and IOV pair.
 *   2. For each record and type that user whant to load, it ask method get for the payload.
 *   3. Reads the data with payloads extracted from step 2.
 *
 * It is not adviced for user to use this class and direct Pool DB manipulation. One should use
 * DataReader and DataWriter classes.
 *
 * The good point to note is that IOV of all L1 trigger condfiguration is controled bay IOV of L1TriggeKey.
 * If new configuration has to be created - new L1TriggerKey has to be saved/loaded. More then one key can use
 * the same paylaod token. This would just mean that data pointed by this payload token has not changed.
 */
class L1TriggerKey
{
public:
    typedef std::map<std::string, std::string> RecordToKey;

    // Constructors
    L1TriggerKey () {}

    /* Adds new record and type mapping to payload. If such exists, nothing happens */
    void add (const std::string & record, const std::string & type, const std::string & key)
    { m_recordToKey.insert (std::make_pair (record + "@" + type, key)); }

    void setTSCKey( const std::string& tscKey )
    { m_tscKey = tscKey ; }

    /* Gets payload key for record and type. If no such paylaod exists, emtpy string
     * is returned.
     */
    std::string get (const std::string & record, const std::string & type) const
    {
        RecordToKey::const_iterator it = m_recordToKey.find (record + "@" + type);
        if (it == m_recordToKey.end ())
            return std::string ();
        else
            return it->second;
    }

    std::string getTSCKey() const
      { return m_tscKey ; }

    const RecordToKey& recordToKeyMap() const
      { return m_recordToKey ; }

protected:
    /* Mapping from records and types to tokens.
     * I as unvable to make type std::map<std::pair<std::string, std::string>, std::string> persistent
     * so record and type are concatanated with @ sign and resulting string is used as a key.
     */

  // wsun 03/2008: instead of tokens, store the configuration keys instead.
/*     typedef std::map<std::string, std::string> RecordsToToken; */
/*     RecordsToToken recordsToToken; */
    RecordToKey m_recordToKey;


    // wsun 03/2008: add data member for TSC key
    std::string m_tscKey ;
};

#endif

