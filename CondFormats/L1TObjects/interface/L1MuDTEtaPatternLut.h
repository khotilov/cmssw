//-------------------------------------------------
//
/**  \class L1MuDTEtaPatternLut
 *
 *   Look-up table for eta track finder
 *
 *
 *   $Date: 2007/02/27 11:44:00 $
 *   $Revision: 1.2 $
 *
 *   N. Neumeister            CERN EP
 */
//
//--------------------------------------------------
#ifndef L1MUDT_ETAPATTERN_LUT_H
#define L1MUDT_ETAPATTERN_LUT_H

//---------------
// C++ Headers --
//---------------

#include <map>

//----------------------
// Base Class Headers --
//----------------------


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

class L1MuDTEtaPattern;

//              ---------------------
//              -- Class Interface --
//              ---------------------

class L1MuDTEtaPatternLut {

  public:

    typedef std::map<int, L1MuDTEtaPattern*, std::less<int> > LUT;
    typedef LUT::const_iterator ETFLut_iter;
    typedef LUT::iterator       ETFLut_Iter;
    
    /// constructor
    L1MuDTEtaPatternLut();

    /// destructor
    virtual ~L1MuDTEtaPatternLut();

    /// reset pattern look-up table
    void reset();
    
    /// load pattern look-up table
    int load();

    /// print pattern look-up table
    void print() const;

    /// get pattern with a given ID
    L1MuDTEtaPattern* getPattern(int id) const;
    
    /// return number of entries in the LUT
    inline int size() const { return m_lut.size(); }

    /// return iterator which points to the first entry of the LUT
    inline ETFLut_iter begin() const { return m_lut.begin(); }

    /// return iterator which points to the one-past-last entry of the LUT
    inline ETFLut_iter end() const { return m_lut.end(); }

  private:

    LUT m_lut;
    
};

#endif
