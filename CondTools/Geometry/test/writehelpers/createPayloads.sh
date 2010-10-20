#!/bin/sh

if [ $# -ne 1 ]
then
  echo Error: createPayloads.sh requires exactly one argument which is the tag
  exit 1
fi
mytag=$1
echo ${mytag}

# Set the tag in all the scripts and the metadata text files
sed -i {s/TagXX/${mytag}/g} *.py
sed -i {s/TagXX/${mytag}/g} *.txt
sed -i {s/TagXX/${mytag}/g} splitDatabase.sh

# First read in the little XML files and create the
# large XML file for the Extended scenario.
# Input cff                    Output file
# GeometryExtended_cff         geSingleBigFile.xml
cmsRun geometryxmlwriter.py

# Now convert the content of the large XML file into
# a "blob" and write it to the database.
# Also reads in the little XML files again and fills
# the DDCompactView. From the DDCompactView the
# reco parts of the database are also filled.
cmsRun geometrywriter.py

# Now put the other scenarios into the database.
# Input the many XML files referenced by the cff file and
# output a single big XML file.
# This is repeated several times below.  The sed commands
# serve to give the following sequence of input and output
# files
#
# Input cff                    Output file
# GeometryExtendedGFlash_cff   gegSingleBigFile.xml
# GeometryIdeal_cff            giSingleBigFile.xml
# GeometryExtendedX0Min_cff    gexminSingleBigFile.xml
# GeometryExtendedX0Max_cff    gexmaxSingleBigFile.xml
# GeometryExtendedLiMin_cff    geliminSingleBigFile.xml
# GeometryExtendedLiMax_cff    gelimaxSingleBigFile.xml
#
sed -i '{s/Extended/ExtendedGFlash/g}' geometryxmlwriter.py
sed -i '{s/\/ge/\/geg/g}' geometryxmlwriter.py
cmsRun geometryxmlwriter.py
sed -i '{s/ExtendedGFlash/Ideal/g}' geometryxmlwriter.py
sed -i '{s/\/geg/\/gi/g}' geometryxmlwriter.py
cmsRun geometryxmlwriter.py
sed -i '{s/Ideal/ExtendedX0Min/g}' geometryxmlwriter.py
sed -i '{s/\/gi/\/gexmin/g}' geometryxmlwriter.py
cmsRun geometryxmlwriter.py
sed -i '{s/X0Min/X0Max/g}' geometryxmlwriter.py
sed -i '{s/\/gexmin/\/gexmax/g}' geometryxmlwriter.py
cmsRun geometryxmlwriter.py
sed -i '{s/X0Max/LiMin/g}' geometryxmlwriter.py
sed -i '{s/\/gexmax/\/gelimin/g}' geometryxmlwriter.py
cmsRun geometryxmlwriter.py
sed -i '{s/LiMin/LiMax/g}' geometryxmlwriter.py
sed -i '{s/\/gelimin/\/gelimax/g}' geometryxmlwriter.py
cmsRun geometryxmlwriter.py

# Read the one big XML file and output a record to the
# database with the an identifying tag
# This is repeated several times below.  The sed commands
# serve to give the following sequence of input file and output
# tag
#
# Input file                Output tag
# gegSingleBigFile.xml      XMLFILE_Geometry_${mytag}_ExtendedGFlash_mc
# giSingleBigFile.xml       XMLFILE_Geometry_${mytag}_Ideal_mc
# gexminSingleBigFile.xml   XMLFILE_Geometry_${mytag}_ExtendedX0Min_mc
# gexmaxSingleBigFile.xml   XMLFILE_Geometry_${mytag}_ExtendedX0Max_mc
# geliminSingleBigFile.xml  XMLFILE_Geometry_${mytag}_ExtendedLiMin_mc
# gelimaxSingleBigFile.xml  XMLFILE_Geometry_${mytag}_ExtendedLiMax_mc
#
sed -i '{s/Extended/ExtendedGFlash/g}' xmlgeometrywriter.py
sed -i '{s/\/ge/\/geg/g}' xmlgeometrywriter.py
cmsRun xmlgeometrywriter.py
sed -i '{s/ExtendedGFlash/Ideal/g}' xmlgeometrywriter.py
sed -i '{s/\/geg/\/gi/g}' xmlgeometrywriter.py
cmsRun xmlgeometrywriter.py
sed -i '{s/Ideal/ExtendedX0Min/g}' xmlgeometrywriter.py
sed -i '{s/\/gi/\/gexmin/g}' xmlgeometrywriter.py
cmsRun xmlgeometrywriter.py
sed -i '{s/X0Min/X0Max/g}' xmlgeometrywriter.py
sed -i '{s/\/gexmin/\/gexmax/g}' xmlgeometrywriter.py
cmsRun xmlgeometrywriter.py
sed -i '{s/X0Max/LiMin/g}' xmlgeometrywriter.py
sed -i '{s/\/gexmax/\/gelimin/g}' xmlgeometrywriter.py
cmsRun xmlgeometrywriter.py
sed -i '{s/LiMin/LiMax/g}' xmlgeometrywriter.py
sed -i '{s/\/gelimin/\/gelimax/g}' xmlgeometrywriter.py
cmsRun xmlgeometrywriter.py

# All the database objects were written into one database
# (myfile.db) in the steps above.  Extract the different
# pieces into separate database files.  These are the payloads
# that get uploaded to the dropbox.  There is one for each tag
./splitDatabase.sh
