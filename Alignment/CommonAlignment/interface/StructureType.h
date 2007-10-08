#ifndef Alignment_CommonAlignment_StructureType_h
#define Alignment_CommonAlignment_StructureType_h

/** \enum StructureType
 *
 *  Enumerate the types of structure an alignable can be.
 *
 *  Basically list the levels in the detector's hierarchy.
 *
 *  $Date: 2007/04/09 00:40:21 $
 *  $Revision: 1.7 $
 *  \author Chung Khim Lae
 */

namespace align
{
  enum StructureType 
  { 
    invalid = 0,
    AlignableDetUnit,
    AlignableDet,

    // Barrel Pixel
    TPBModule,
    TPBLadder,
    TPBLayer,
    TPBHalfBarrel,
    TPBBarrel,

    // Forward Pixel
    TPEModule,
    TPEPanel,
    TPEBlade,
    TPEHalfDisk,
    TPEHalfCylinder,
    TPEEndcap,

    // Tracker Inner Barrel
    TIBModule,
    TIBString,
    TIBSurface,
    TIBHalfShell,
    TIBLayer,
    TIBHalfBarrel,
    TIBBarrel,

    // Tracker Inner Disks
    TIDModule,
    TIDSide,
    TIDRing,
    TIDDisk,
    TIDEndcap,

    // Tracker Outer Barrel
    TOBModule,
    TOBRod,
    TOBLayer,
    TOBHalfBarrel,
    TOBBarrel,

    // Tracker Endcaps
    TECModule,
    TECRing,
    TECPetal,
    TECSide,
    TECDisk,
    TECEndcap,

    Pixel,
    Strip,
    Tracker,

    // Muon Detector, not touching these now
    AlignableDTBarrel = 40,
    AlignableDTWheel,
    AlignableDTStation,
    AlignableDTChamber,
    AlignableDTSuperLayer,
    AlignableDTLayer,
    AlignableCSCEndcap,
    AlignableCSCStation,
    AlignableCSCChamber,
    AlignableCSCLayer,
    AlignableMuon,

    Detector
  };
}

#endif
